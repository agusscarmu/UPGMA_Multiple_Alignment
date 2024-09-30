import numpy as np
from Bio import SeqIO
import matplotlib.pyplot as plt
from Bio.Align import substitution_matrices
from Bio import Phylo
from io import StringIO
import pandas as pd
import argparse
import time
import threading
import warnings
import platform
import os

# Suprimir el UserWarning específico
warnings.filterwarnings("ignore", category=UserWarning, message="FigureCanvasAgg is non-interactive, and thus cannot be shown")


def read_fasta(filename):
    with open(filename,'r') as handle:
        r_seq = list(SeqIO.parse(handle, 'fasta'))
        return r_seq[0].seq, r_seq[1].seq
    
class Alignment: 
    def __init__(self, seq1, seq2, k, blosum, type_of_sequences, threshold=0):
        """
        Genera la mejor alineación entre dos secuencias utilizando un dotplot filtrado.
        
        Args:
            seq1 (str): Primera secuencia a alinear.
            seq2 (str): Segunda secuencia a alinear.
            k (int): Tamaño de las submatrices cuadradas utilizadas en el filtrado del dotplot.
            t (str): Tipo de secuencia ('PROTEIN' o 'DNA').
        
        Returns:
            tuple: Una tupla que contiene la alineación traducida, el dotplot original, y el dotplot filtrado.
        
        Raises:
            ValueError: Si el tamaño `k` es mayor que las dimensiones de la matriz `dotplot`.
        """
        if len(seq1[0])<2 and len(seq2[0])<2:
            raise ValueError(f"El formato de la secuencia debe ser: [[sequence],'name of sequence']")
        if len(seq1[0]) <= len(seq2[0]):
            self.seq1 = seq1[0]
            self.seq2 = seq2[0]
        else:
            self.seq1 = seq2[0]
            self.seq2 = seq1[0]
        self.type_of_sequences = type_of_sequences.upper()
        self.blosum = blosum
        self.seq1_label = seq1[1]
        self.seq2_label = seq2[1]
        if self.type_of_sequences == 'DNA':
            self.threshold = (k/2)
        elif threshold==0:
            self.threshold = self.__obtain_threshold()          
        if k > len(seq1[0]) or k > len(seq2[0]):
            raise ValueError(f"El tamaño de {k} debe ser menor al tamaño de la secuencia.")
        self.k = k

        self.dotplot = self.__make_dotplot(self.seq1, self.seq2)
        self.f_dotplot = self.__filtered_dotplot(k)
        self.best_alignment = self.__obtain_best_alignment()
        self.traduced_alignment = self.__traduce_alignment()
        self.score = self.__get_score()

    def __get_score(self):
        score = 0.0
        for coord in self.best_alignment:
            score += self.dotplot[coord[0]][coord[1]]
        return score
    
    def __obtain_threshold(self):
        t_dotplot1, t_dotplot2 = self.__make_dotplot(self.seq1, self.seq1), self.__make_dotplot(self.seq2, self.seq2)
        i=0
        j=0
        threshold1, threshold2 = 0, 0
    
        while i < len(t_dotplot1) and j < len(t_dotplot1[0]):
            threshold1 += t_dotplot1[i, j]
            i+=1
            j+=1
        i=0
        j=0
        while i < len(t_dotplot2) and j < len(t_dotplot2[0]):
            threshold2 += t_dotplot2[i, j]
            i+=1
            j+=1
        threshold = ((threshold1/len(self.seq1)) + (threshold2/len(self.seq2))/2)
    
        return threshold

    def __traduce_alignment(self):
        traduced_seq = ''
        for i in range(len(self.best_alignment)):
            if self.seq1[self.best_alignment[i][0]] == self.seq2[self.best_alignment[i][1]]:
                traduced_seq += f"{self.seq1[self.best_alignment[i][0]]}{i}"
            else:
                # traduced_seq += f"{seq1[alignment_list[i-1][0]]}"
                # Si se prefiere un guion en vez de repetir la ultima base:
                traduced_seq += "-"
        return traduced_seq

    def __obtain_best_alignment(self):
        complete = False
        alignment_list = []
        i = 0
        j = 0
    
        while not complete:
            # Verificar si hemos llegado al final de alguna secuencia
            if i >= len(self.f_dotplot) or j >= len(self.f_dotplot[0]):
                complete = True
            else:
                alignment_list.append([i, j])  # Agregar la posición actual al alineamiento
    
                # Verificar si podemos movernos en diagonal
                if i + 1 < len(self.f_dotplot) and j + 1 < len(self.f_dotplot[0]) and self.f_dotplot[i + 1, j + 1] == 1:
                    i += 1
                    j += 1
                # Verificar si podemos movernos horizontalmente (gap en seq1)
                elif j + 1 < len(self.f_dotplot[0]) and self.f_dotplot[i, j + 1] == 1:
                    j += 1
                # Verificar si podemos movernos verticalmente (gap en seq2)
                elif i + 1 < len(self.f_dotplot) and self.f_dotplot[i + 1, j] == 1:
                    i += 1
                else:
                    if i<len(self.f_dotplot) and j<len(self.f_dotplot[0]):
                        j+=1
                    else:
                        complete = True
    
        return alignment_list



    def __has_potential(self, subdotplot):
        score = 0
        for i in range(len(subdotplot)):
            score += subdotplot[i,i]
    
        if score > self.threshold:
            return True
        return False
    
    def __filtered_dotplot(self,k):
        """
        Filtra una matriz de dotplot resaltando las regiones diagonales de tamaño k con mayor puntuación.
        
        La función evalúa las submatrices de tamaño `k x k` dentro de la matriz `dotplot`, buscando aquellas que tienen un número mayoritario de unos en su diagonal principal. Si la submatriz cumple con el criterio, se añade al dotplot filtrado.
        
        Returns:
            numpy.ndarray: Una nueva matriz del mismo tamaño que `dotplot`, donde las regiones que cumplen con 
            el criterio están resaltadas con unos en la diagonal.
        
        """
        rows, cols = len(self.dotplot), len(self.dotplot[0])
        f_dotplot = np.zeros((rows, cols))
    
        for i in range(rows-k+1):
            for j in range(cols-k+1):
                sub_dotplot = self.dotplot[i:i+k, j:j+k]
                if self.__has_potential(sub_dotplot):
                    if self.type_of_sequences == 'DNA':
                        f_dotplot[i:i+k, j:j+k] = sub_dotplot * np.eye(k)
                    elif self.type_of_sequences == 'PROTEIN':
                        f_dotplot[i,j]=1
    
        return f_dotplot

    def show_dotplot(self, filename=None, filtered=False, show_values=False, show_alignment=False):
        
        if filtered:
            dtplt = self.f_dotplot
        else:
            dtplt = self.dotplot
            
        plt.imshow(dtplt, cmap="gray_r", interpolation='nearest')
        if not filtered:
            k_message = "SIN FILTRAR"     
        else:
            k_message = f"k={self.k}"
    
        if show_values:
            for i in range(len(dtplt)):
                for j in range(len(dtplt[0])):
                    value = dtplt[i][j]
                    plt.text(j, i, f'{value:.1f}', ha='center', va='center', color='red', fontsize=8)

        num_rows = len(dtplt)
        num_cols = len(dtplt[0])
        if num_rows > num_cols:
            point_size = num_rows/num_cols
        else: 
            point_size = num_cols/num_rows
        # Mostrar el alineamiento en rojo (opcional)
        if show_alignment and hasattr(self, 'best_alignment'):
            for alg in self.best_alignment:
                # Limitar los puntos a las dimensiones del dotplot
                if 0 <= alg[0] < num_rows and 0 <= alg[1] < num_cols:
                    plt.scatter(alg[1], alg[0], color='red', s=point_size)  # Trazar los puntos del alineamiento
                
        plt.title(f'Dotplot: {len(dtplt)} elementos vs {len(dtplt[0])} elementos | {k_message}')
        plt.xlabel(self.seq2_label)
        plt.ylabel(self.seq1_label)
        plt.show()
        if filename:
            plt.savefig(filename)
    
    def __prepare_seq_to_show_alignment(self, seq, position, len_of_seq2, len_score):
        prepared_seq = ""
        i=0
        while i < len_of_seq2:
            if i == position:
                j=0
                while j < len_score:
                    prepared_seq += seq[j]
                    j+=1
                i=i+j
            else:
                prepared_seq += "-"
            i+=1
        return prepared_seq
    
    def seq_string_equals(self, seq1, seq2):
        if self.seq1_label == seq1 and self.seq2_label == seq2:
            return True
        elif self.seq1_label == seq2 and self.seq2_label == seq1:
            return True
        else:
            return False
        
    def show_alignments(self):
        match_line = ''
        match_score = 0
        best_score = 0
        position = 0
        
        for i in range(len(self.best_alignment)):
            base_seq1 = self.seq1[self.best_alignment[i][0]]
            base_seq2 = self.seq2[self.best_alignment[i][1]]
    
            # Verificar si son iguales
            if base_seq1 == base_seq2:
                match_line += '|'
                match_score += 1
            else:
                match_line += '.'
                if match_score > best_score:
                    best_score = match_score
                    position = i-best_score
                match_score = 0

        print(self.__prepare_seq_to_show_alignment(self.seq1, position, len(self.seq2), best_score))
        print(match_line)
        print(self.seq2)
    
    def __make_dotplot(self, s1, s2):
        dotplot = np.zeros((len(s1),len(s2)))
    
        for i in range(len(s1)):
            for j in range(len(s2)):
                if self.type_of_sequences=='DNA':
                    if s1[i]==s2[j]:
                        dotplot[i,j] = 1
                else:
                    dotplot[i,j] = substitution_matrices.load(self.blosum)[s1[i],s2[j]]
    
        return dotplot    
    

class AlignmentNode:

    def __init__(self, seq1, seq2):
        self.seq1 = seq1
        self.seq2 = seq2
        self.__list_of_traduced_alignments = []

    def __obtain_list_seqs(self, seq, label):
        k = len(label)
        if isinstance(seq, str):
            l = []
            i = 0
            while i < len(seq):
                if seq[i:i+k] == label:  # Encuentra 'seq'
                    j = i + k
                    # Encuentra el número después de 'seq'
                    while j < len(seq) and seq[j].isdigit():
                        j += 1
                    l.append(seq[i:j])  # Captura 'seq' seguido del número
                    i = j 
                else:
                    i += 1  # Avanza si no es 'seq'
            return l
        elif isinstance(seq, AlignmentNode):
            return seq.__obtain_list_seqs(seq.seq1,label)+(seq.__obtain_list_seqs(seq.seq2,label))
    
    # Función recursiva para traer los alineamientos
    def __recursive_alignment(self, seq, d_matrix, alignment):
        # Primero recorre de manera recursiva (excepto si la rama es str) trayendo la lista de mejores secuencias
        if isinstance(seq.seq1, str):
            seq_to_align_1 = [seq.seq1]
        else:
            seq_to_align_1 = self.__recursive_alignment(seq.seq1, d_matrix, alignment)
        if isinstance(seq.seq2, str):
            seq_to_align_2 = [seq.seq2]
        else:
            seq_to_align_2 = self.__recursive_alignment(seq.seq2, d_matrix, alignment)
        
        # Controla las secuencias traidas (todas contra todas) y se fija cual es la de menor distancia (mejor score)
        best_score = -np.inf
        best_seqs = []
        for i in range(len(seq_to_align_1)):
            for j in range(len(seq_to_align_2)):
                if d_matrix.loc[seq_to_align_1[i]][seq_to_align_2[j]] > best_score:
                    best_score = d_matrix.loc[seq_to_align_1[i]][seq_to_align_2[j]]
                    best_seqs = [seq_to_align_1[i], seq_to_align_2[j]]
        
        # Recorre la lista de alineamientos y busca el alineamiento que coincide con el mejor score obtenido previamente
        for algn in alignment:
            if algn.seq_string_equals(best_seqs[0], best_seqs[1]):
                self.__list_of_traduced_alignments.append([algn.traduced_alignment, best_seqs])
                break
                
        return best_seqs
            
    def alignment(self, d_matrix, alignments):
        self.__recursive_alignment(self, d_matrix, alignments)
        return self.__list_of_traduced_alignments.copy()
    
    def contains_seq(self, seq, label):
        list = self.__obtain_list_seqs(self, label)
        list2 = self.__obtain_list_seqs(seq, label)
        for seq1 in list: 
            for seq2 in list2:
                if seq1 == seq2:
                    return True
        return False
            
    
    def __str__(self):
        return f"({self.seq1}, {self.seq2})"


class MultipleAlignment:   
    def __init__(self, sequences, blosum, type_of_sequences="DNA"):
        # Lista de nodos de alineamiento
        self.list_of_nodes = []
        # Alineamientos totales (util para graficar el arbol)
        self.alignments = None
        # Historial de matriz de distancia, permite controlar la iteracion completa
        self.matrix_history = []
        # Nombre de las secuencias
        self.label_of_sequences = ""
        # Matriz de distancia inicial
        self.distance_matrix = None
        # Alineamientos iniciales
        self.sequences_alignments = []
        self.blosum = blosum
        # Secuencias a alinear
        self.sequences = sequences
        type_of_sequences = type_of_sequences.upper()
        if type_of_sequences == "DNA" or type_of_sequences == "PROTEIN":     
                self.type_of_sequences = type_of_sequences
                if self.type_of_sequences == "DNA":
                    self.blosum = None
        else:
            raise ValueError("type_of_sequences debe ser DNA o PROTEIN")
        if self.__check_sequences():        
            self.matrix = self.__obtain_matrix()
        else:
            raise ValueError("El formato de la secuencia tiene que ser:  [[sequence],'name of sequence']")

    def __check_sequences(self):
        # Verificar que la lista contenga tuplas o listas con dos elementos
        if not isinstance(self.sequences, list):
            return False
    
        for item in self.sequences:
            if not isinstance(item, (tuple, list)) or len(item) != 2:
                return False
    
            sequence, name = item
            
            if self.type_of_sequences == "DNA":
                # Verificar que la secuencia sea una cadena que solo contenga A, T, C, G
                if not isinstance(sequence, str) or not all(base in 'ATCG' for base in sequence.upper()):
                    raise ValueError(f"La secuencia es: '{self.type_of_sequences}' pero no contiene 'ATCG'")
                
            # Verificar que el nombre sea una cadena
            if not isinstance(name, str):
                raise ValueError(f"El nombre de la secuencia tiene que ser un String")
    
        return True

        
    # Dibuja el dendograma
    def draw_tree(self, filename=None):  
        tree = Phylo.read(StringIO(f"{self.alignments}"), "newick")
        tree.ladderize()
        if filename:
            Phylo.draw(tree, do_show=False)
            plt.savefig(filename)
        Phylo.draw(tree)
    
    # Muestra todos los alineamientos respetando el dendograma
    def show_alignment(self):
        for alignment in self.alignments.alignment(self.distance_matrix, self.sequences_alignments):
            print(f"Secuencias: {alignment[1][0], alignment[1][1]}\n"
                  f"Alineamiento: {alignment[0]}\n"
                  f"{'_'*30}")
    
    def obtain_alignment(self, seq1, seq2):
        for alignment in self.sequences_alignments:
            if alignment.seq_string_equals(seq1, seq2):
                return alignment
        
        return None
    
    # Encuentra las secuencias con score mas cercano
    def __closest_neighbor(self,matrix):
        pair = []
        best_distance = -np.inf
        for i in range(len(matrix.index)):
            for j in range(i+1,len(matrix.columns)):
                distance = matrix.loc[matrix.index[i], matrix.columns[j]]
                if distance > best_distance:
                    pair=[matrix.index[i], matrix.columns[j]]
                    best_distance = matrix.loc[matrix.index[i], matrix.columns[j]]
        return pair[0], pair[1]
    
    # Revisa si existe el nodo, de ser asi lo devuelve
    def __check_if_exists_node(self, seq):
        r = None
        for node in self.list_of_nodes:
            if node.contains_seq(seq, self.label_of_sequences):
                r = node
            
        if r is None:
            return seq
        else:
            self.list_of_nodes.remove(r)
            return r
    
    # Se edita la matriz eliminando las columnas y las filas de los nodos conectados, y creando una nueva con el nodo nuevo.
    def __edit_matrix(self,matrix, row, column):
        seq1 = row
        seq2 = column
        
        df_copy = matrix.copy()
        matrix = matrix.drop(index=seq1, columns=seq2)
        matrix = matrix.drop(index=seq2, columns=seq1)
        
        align_node = AlignmentNode(self.__check_if_exists_node(seq1), self.__check_if_exists_node(seq2))
        self.list_of_nodes.append(align_node)
        # Insertar una columna llena de ceros en la posición dinámica
        
        matrix.insert(0, f'{align_node}', [0] * len(matrix))
        
        # Crear una fila llena de ceros para insertar con un nombre específico
        new_row = pd.DataFrame([[0] * len(matrix.columns)], columns=matrix.columns, index=[f'{align_node}'])
        
        # Insertar la fila en la posición dinámica
        matrix = pd.concat([new_row, matrix.iloc[0:]])

        matrix = matrix.astype(float)
        
        node1_str = str(align_node.seq1) 
        node2_str = str(align_node.seq2)
        
        # Reasigna valores a la matriz (en las filas y columnas nuevas hace la formula mencionada en el markdown de arriba)
        for i in range(len(matrix.index)):
            for j in range(i+1,len(matrix.columns)): 
                if matrix.index[i] == f"{align_node}": 
                    matrix.loc[matrix.index[i],matrix.columns[j]] = (df_copy.loc[node1_str,matrix.columns[j]]  + df_copy.loc[matrix.columns[j],node2_str])/2
                elif matrix.columns[j] == f"{align_node}":
                    matrix.loc[matrix.index[i],matrix.columns[j]] = (df_copy.loc[node1_str,matrix.index[i]]  + df_copy.loc[matrix.index[i],node2_str])/2
                else:
                    matrix.loc[matrix.index[i],matrix.columns[j]] = df_copy.loc[matrix.index[i],matrix.columns[j]]        
        return matrix, align_node
    
    
    # Obtiene la matriz de distancia
    def __d_matrix(self):
        n = len(self.sequences)
        matrix = np.zeros((n, n))
        for i in range(n):
            for j in range(i+1, n):
                new_align = Alignment(self.sequences[i], self.sequences[j], 3, self.blosum, self.type_of_sequences, )
                self.sequences_alignments.append(new_align)
                matrix[i, j] = new_align.score
        df_distance_matrix = pd.DataFrame(matrix)
        self.label_of_sequences = self.sequences[0][1][:-1]
        df_distance_matrix.columns = [x[1] for x in self.sequences]
        df_distance_matrix.index = [x[1] for x in self.sequences]

        return df_distance_matrix
    
    # Función principal de la clase, mediante esta va obteniendo la matriz y generando los nodos
    def __obtain_matrix(self):
        # Obtiene la matriz de distancia de las secuencias (compara todas con todas)
        matrix = self.__d_matrix()
        self.distance_matrix = matrix
        while len(matrix)>1:
            self.matrix_history.append(matrix)
            row, column = self.__closest_neighbor(matrix)
            matrix, align_node = self.__edit_matrix(matrix, row, column)
            self.alignments = align_node

def read_fasta_to_list(file_path):
    sequences = []
    with open(file_path, 'r') as file:
        seq_name = None
        seq_data = []

        for line in file:
            line = line.strip()
            if line.startswith('>'):  
                if seq_name:  
                    sequences.append([ ''.join(seq_data), seq_name ])
                seq_name = line[1:]  
                seq_data = []  
            else:
                seq_data.append(line)  

        # Guardar la última secuencia después de terminar el bucle
        if seq_name:
            sequences.append([ ''.join(seq_data), seq_name ])

    return sequences

stop_event = threading.Event()


def loading():
    while not stop_event.is_set():
        for indicador in ["|", "/", "-", "\\"]:
            if stop_event.is_set():
                break
            print(f"\rAligning sequences... {indicador}", end="")
            time.sleep(0.2)

def clear_console():
    if platform.system() == "Windows":
        os.system("cls")
    else:
        os.system("clear")

def start_multiple_alignment(fasta, type_of_sequences, blosum):
    thread_loading = threading.Thread(target=loading)
    thread_loading.daemon = True
    thread_loading.start()

    sequences = read_fasta_to_list(fasta)
    multiple_alignment = MultipleAlignment(sequences, blosum, type_of_sequences)

    stop_event.set()
    thread_loading.join()

    option = 0


    while option != 4:
            clear_console()
            print("\nMULTIPLE SEQUENCE ALIGNMENT\n------------------------------------------------------")
            print("\n1. Show dendogram")
            print("2. Show alignments")
            print("3. Get alignment between two sequences")
            print("4. Exit")
            while True:
                try:
                    option = int(input("Enter an option: "))
                    break
                except ValueError:
                    print("Please enter a valid number.")

            if option == 1:
                file_name = input("Enter the file name to save the dendogram: ")
                directory = "graphics/dendograms/"
                
                # Verificar si el directorio existe y crearlo si no es así
                if not os.path.exists(directory):
                    os.makedirs(directory)
                
                file_name = directory + file_name + ".png"
                multiple_alignment.draw_tree(filename=file_name)  
                print(f"The dendogram has been saved in {file_name}")
                time.sleep(2)  
            elif option == 2:
                multiple_alignment.show_alignment()
                input("Press any key to continue...")
            elif option == 3:
                seq1 = input("Enter the name of the first sequence: ")
                seq2 = input("Enter the name of the second sequence: ")
                alignment = multiple_alignment.obtain_alignment(seq1, seq2)
                if alignment is None:
                    print("No alignment was found between the entered sequences.")
                    time.sleep(2)
                else:
                    option_3 = 0
                    while option_3 != 5:
                        clear_console()
                        print(f"\nALIGNMENT OF {seq1.upper()} AND {seq2.upper()}\n------------------------------------------------------")
                        print("\n1. Show dotplot")
                        print("2. Show filtered dotplot")
                        print("3. Show alignment")
                        print("4. Get Score")
                        print("5. Back to multiple alignment options")

                        while True:
                            try:
                                option_3 = int(input("Enter an option: "))
                                break
                            except ValueError:
                                print("Please enter a valid number.")

                        if option_3 == 1 or option_3 == 2:
                            file_name = input("Enter the file name to save the dotplot: ")
                            directory = "graphics/dotplots/"

                            # Verificar si el directorio existe y crearlo si no es así
                            if not os.path.exists(directory):
                                os.makedirs(directory)

                            file_name = directory+file_name + ".png"
                            alignment.show_dotplot(filename=file_name) if option_3 == 1 else alignment.show_dotplot(filename=file_name, filtered=True, show_alignment=True)
                            print(f"The dotplot has been saved in {file_name}")
                            time.sleep(2) 
                        elif option_3 == 3:
                            print(alignment.traduced_alignment)
                            input("Press any key to continue...")
                        elif option_3 == 4:
                            print(f"Score: {alignment.score}")
                            input("Press any key to continue...")

    

    

def main():
    # Crear el parser
    parser = argparse.ArgumentParser(description='Script para alineamiento múltiple de secuencias.')

    # Argumento posicional: archivo FASTA
    parser.add_argument('fasta_file', type=str, help='Archivo de secuencias en formato FASTA.')

    parser.add_argument('type_of_sequences', type=str, choices=['DNA', 'PROTEIN'],help='Tipo de secuencias (DNA o PROTEIN).')

    # Argumento opcional: matriz BLOSUM
    parser.add_argument('--blosum', type=str, default='BLOSUM62', 
                        help='Especificar la matriz BLOSUM a usar (default: BLOSUM62).')

    # Parsear los argumentos
    args = parser.parse_args()

    start_multiple_alignment(args.fasta_file, args.type_of_sequences, args.blosum)


if __name__ == "__main__":
    main()