import csv
import os
import ast
import numpy as np
import seaborn as sns
import scipy
import matplotlib.pyplot as plt
import sys
from kuramoto import Kuramoto, plot_activity, plot_phase_coherence
from sklearn.decomposition import PCA
from tqdm import tqdm
from iced import normalization
from tqdm.contrib.concurrent import process_map
import concurrent
import time
import multiprocessing as mp
from multiprocessing import Pool
from itertools import repeat
    
def submatrix_reader(matrix,submatrix_shape,stride = 1,offset = (0,0),pattern = "diagonal"):
    #Logic should be working only for symmetric submatrices, not properly generalized to non square matricess
    submatrix_no = min(matrix.shape) - max(submatrix_shape) + 1 - max(offset)
    submatrices=[]
    for i in range(0,submatrix_no,stride):
        submatrix = matrix[i+offset[0]:i+submatrix_shape[0]+offset[0],i+offset[1]:i+submatrix_shape[1] + offset[1]]
        submatrices.append(submatrix)
    return submatrices

def universal_order_parameter(phase_time_series_matrix,adj_mat):
        '''
        Compute universal order parameter R_t - mean length of resultant vector
        '''
        suma=0
        angles_vec = phase_time_series_matrix
        for i in range(0,adj_mat.shape[0]):
            for j in range(0,adj_mat.shape[1]):
                #if i>j:
                    #continue
                i_vec = angles_vec[i][~np.isnan(angles_vec[i])]
                j_vec = angles_vec[j][~np.isnan(angles_vec[j])]
                if len(i_vec) > 1000:
                     suma += adj_mat[i,j]*sum(np.cos(i_vec[-1000:]-j_vec[-1000:]))/len(i_vec[-1000:])
                     #print(suma,sum(np.cos(i_vec-j_vec))/1000*adj_mat[i,j],adj_mat[i,j],i_vec-j_vec)
                else:
                    suma += adj_mat[i,j]*sum(np.cos(i_vec-j_vec))/len(i_vec)
        denom = np.sum(adj_mat)
        #print(suma,denom)
        return abs(suma / denom)

def list_writer(filename,variable_name):
    file = open('{}'.format(filename), 'w')
    for element in variable_name:
            file.write(str(element) +'\t')
    file.close()

def matrix_reader(filename,delete = False):
    file = open('{}'.format(filename), 'r')
    content = file.readlines()
    #print(len(content[1]))
    file.close()

    tsv_reader = csv.reader(content, delimiter='\t')
    counter = 0
    row_list = []
    for row in tsv_reader :
            if delete == True:
                del row[len(row)-1]
            for i in range(0,len(row)):
                 row[i] = float(row[i])
            counter += 1
            row_list.append(row)
    return np.array(row_list)

def matrix_writer(filename,variable_name):
    file = open('{}'.format(filename), 'w')

    for row in variable_name:
        for element in row:
            file.write(str(element) +'\t')
        file.write('\n')
    file.close()      

def z_normalization(matrix,iteration_no = 10):
    #Numpy arrays are mutable so we first make a copy to preserve the original array and then use mutability to our advantage
    matrix = matrix.copy()
    for _ in range(0,iteration_no):
        for i in range(0,matrix.shape[0]):
            matrix[i,:] = (matrix[i,:]- np.mean(matrix[i,:]))/(np.std(matrix[i,:])+1e-48)
        for j in range(0,matrix.shape[1]):
            matrix[:,j] = (matrix[:,j]- np.mean(matrix[:,j]))/(np.std(matrix[:,j])+1e-48)
    return matrix


def k_r_calculator_parallel(HiC_submatrix_entry,no_nodes):
        #k_50_storer_1 = []
        #k_50_storer_2 =[]
        coupling_vals = np.linspace(0, 0.6, 200)
        n_nodes=no_nodes
        #graph = HiC_50[val]
        graph = HiC_submatrix_entry
        adjmat=graph
        runs = []
        for coupling in tqdm(coupling_vals):
            model = Kuramoto(coupling=coupling, dt=0.1, T=500, n_nodes=n_nodes)
            model.natfreqs = np.random.normal(1, 0.1, size=n_nodes)  # reset natural frequencies
            act_mat = model.run(adj_mat=adjmat)
            runs.append(act_mat)

        # Check that natural frequencies are correct (we need them for prediction of Kc)
        plt.figure()
        plt.hist(model.natfreqs)
        plt.xlabel('natural frequency')
        plt.ylabel('count')
        start_time = time.perf_counter()
        # Plot all time series for all coupling values (color coded)
        runs_array = np.array(runs)
        end_time = time.perf_counter()
        elapsed_time = end_time - start_time
        print(f"Execution time: {elapsed_time:.6f} seconds")
        #plt.show()
        plt.close()
        r_k_storer_1 = []
        r_k_storer_2 = []
        plt.figure()
        for i, coupling in tqdm(enumerate(coupling_vals)):
            r_mean_kura = np.mean([model.phase_coherence(vec)
                          for vec in runs_array[i, :, -1000:].T]) # mean over last 1000 steps
            r_mean_uni = universal_order_parameter(runs_array[i],adjmat)
            r_k_storer_1.append((coupling,r_mean_kura))
            r_k_storer_2.append((coupling,r_mean_uni))
            plt.scatter(coupling, r_mean_kura, c='steelblue', s=20, alpha=0.7,label = "Kuramoto")
            plt.scatter(coupling, r_mean_uni, c='red', s=20, alpha=0.7,label = "Universal")

        # Predicted Kc – analytical result (from paper)
        #Kc = np.sqrt(8 / np.pi) * np.std(model.natfreqs) # analytical result (from paper)
        #plt.vlines(Kc, 0, 1, linestyles='--', color='orange', label='analytical prediction')
        #plt.legend()
        #plt.grid(linestyle='--', alpha=0.8)
        #plt.title('Critical Coupling  1000 TAD {} Universal r '.format(val))
        #plt.ylabel('order parameter (R)')
        #plt.xlabel('coupling (K)')
        #sns.despine()
        #plt.show()
        #if val <= 25:
        #plt.savefig('C:/Users/HP/Downloads/Siddharth/Thesis_Figures/Weighted_Network/TAD_50_human_chr11/Normalized/ICE_normalized/Combined/TAD_50_human_chr11_{}_normalized_combined'.format(val))
        plt.close()
        #Storing k and r values for every coupling value
        r_storer_1 = []
        k_storer_1 = []
        for i,k_r in enumerate(r_k_storer_1):
            r_storer_1.append(k_r[1])
            k_storer_1.append(k_r[0])
        r_storer_2 = []
        k_storer_2 = []
        for i,k_r in enumerate(r_k_storer_2):
            r_storer_2.append(k_r[1])
            k_storer_2.append(k_r[0])
        return (r_storer_1,k_storer_1,r_storer_2,k_storer_2,np.sum(HiC_submatrix_entry))
#PCA for compartments
def PCA_application(matrix,no_components = 3):
        pca = PCA(n_components = no_components)
        PCA_output = pca.fit(matrix)
        return PCA_output.components_

#Saving the column sum
def Col_Row_sum_calc(matrix):
        Column_sum_list = []
        Row_sum_list = []
        for i in range(0,matrix.shape[1]):
            Column_sum_list.append(np.sum(matrix[:,i]))
        for j in range(0,matrix.shape[0]):
            Row_sum_list.append(np.sum(matrix[i,:]))
        return (Column_sum_list,Row_sum_list)

#Obtaining insulation score
def insulation_score_calc(HiC_matrix,node_no=10,stride = 1,offset = (0,10)):
        insulated_sub =  submatrix_reader(HiC_matrix,(node_no,node_no),stride,offset)
        insulation_score = list(map(np.sum,insulated_sub))
        return insulation_score

#Obtaining local distal interaction ratio:
def local_distal_interaction_ratio(HiC_matrix):
        Distal_Interaction_List = [0]*HiC_matrix.shape[0]
        Local_Interaction_List = [0]*HiC_matrix.shape[0]

        for i in range(0,HiC_matrix.shape[0]):
            for j in range(0,HiC_matrix.shape[1]):
                if abs(j-i) < 3000000/Resolution :
                    Local_Interaction_List[i]+=HiC_data_to_work[i,j]
                else:
                    Distal_Interaction_List[i]+=HiC_data_to_work[i,j]
        log_2_LDR_list = [0]*HiC_data_to_work.shape[0]
        for i in range(0,HiC_data_to_work.shape[0]):
            if Distal_Interaction_List[i] == 0 :
                Distal_Interaction_List[i] += 0.000001
            log_2_LDR_list[i] = np.log2(Local_Interaction_List[i]/(Distal_Interaction_List[i]))
        return log_2_LDR_list

def multi_matrix_k_r_parallel_part_1(HiC_data_to_work,image_root,data_root,Matrix_Form,node_no=10,processes_to_use = 2):
        #print("Reached")
        if Matrix_Form == "Raw":
            HiC_data_to_work = HiC_data_to_work
        elif Matrix_Form == "ICE":
            HiC_data_to_work = normalization.ICE_normalization(HiC_data_to_work)
        elif Matrix_Form == "Z":
            HiC_data_to_work = z_normalization(HiC_data_to_work)
        hmap = plt.imshow(HiC_data_to_work, cmap='Reds', interpolation='nearest')
        plt.colorbar(hmap)
        Img_Dir_exist = os.path.exists('{}/{}'.format(image_root,Matrix_Form))
        Data_Dir_exist = os.path.exists('{}/{}'.format(data_root,Matrix_Form))
        if not Img_Dir_exist:
            os.makedirs('{}/{}'.format(image_root,Matrix_Form))
        if not Data_Dir_exist:
            os.makedirs('{}/{}'.format(data_root,Matrix_Form))
        plt.savefig('{}/{}/{}_heatmap'.format(image_root,Matrix_Form,Matrix_Form))
        #plt.show()
        plt.close()
        matrix_writer("{}/{}/{}_PCs".format(data_root,Matrix_Form,Matrix_Form),PCA_application(HiC_data_to_work))
        list_writer("{}/{}/{}_Col_sum".format(data_root,Matrix_Form,Matrix_Form),Col_Row_sum_calc(HiC_data_to_work)[0])
        list_writer("{}/{}/{}_Insulation_score".format(data_root,Matrix_Form,Matrix_Form),insulation_score_calc(HiC_data_to_work))
        list_writer("{}/{}/{}_LDR_list".format(data_root,Matrix_Form,Matrix_Form),local_distal_interaction_ratio(HiC_data_to_work))
        HiC_sub = submatrix_reader(HiC_data_to_work,(node_no,node_no))
        return HiC_sub

        #print(HiC_sub[0:10],np.sum(HiC_sub[0]))
        '''
        r_super_storer_1 = []
        k_super_storer_1 = []
        r_super_storer_2 = []
        k_super_storer_2 = []
        submatrix_sum_list = []
        max_ = len(HiC_sub)
        with Pool(processes=processes_to_use) as p, tqdm(total=max_) as pbar:
                for result in p.starmap(k_r_calculator_parallel, zip((HiC_sub),repeat(node_no))):
                    pbar.update()
                    pbar.refresh()
                    r_super_storer_1.append(result[0])
                    k_super_storer_1.append(result[1])
                    r_super_storer_2.append(result[2])
                    k_super_storer_2.append(result[3])
                    submatrix_sum_list.append(result[4])
        matrix_writer("{}/{}/{}_k_values_kura_order".format(data_root,Matrix_Form,Matrix_Form),k_super_storer_1)
        matrix_writer("{}/{}/{}_r_values_kura_order".format(data_root,Matrix_Form,Matrix_Form),r_super_storer_1)
        matrix_writer("{}/{}/{}_k_values_uni_order".format(data_root,Matrix_Form,Matrix_Form),k_super_storer_2)
        matrix_writer("{}/{}/{}_r_values_uni_order".format(data_root,Matrix_Form,Matrix_Form),r_super_storer_2)
        list_writer("{}/{}/{}_submatrix_sum_list".format(data_root,Matrix_Form,Matrix_Form),submatrix_sum_list)
        #p.close()
        #p.join()
        #time.sleep(60)       
        '''
if __name__ == '__main__':
    '''
    if len(sys.argv) > 1 :
        if sys.argv[1] == "system":
            #print("yay")
            HiC_data = matrix_reader(sys.argv[2])
            np.fill_diagonal(HiC_data,0)
            Regions_to_check = ast.literal_eval(sys.argv[3])
            #print(Regions_to_check,type(Regions_to_check),type(Regions_to_check[0]))
            Resolution = sys.argv[4]
            node_no = sys.argv[5]
            data_root = sys.argv[6]
            image_root = sys.argv[7]
    else:
        
            HiC_data = matrix_reader("C:/Users/HP/Downloads/Siddharth/Data/hg38/GM12878/100kb/raw/chr1.mat.txt")
            np.fill_diagonal(HiC_data,0) #Making diagonals of our matrix zero
            Regions_to_check = (143000001,158000001)
            Resolution = 100000
            node_no = 10
            data_root = "C:/Users/HP/Downloads/Siddharth/Data/hg38/GM12878/100kb/Chromosome_1"
            image_root = "C:/Users/HP/Downloads/Siddharth/Thesis_Figures/Weighted_Network/GM12878/100kb/Chromosome_1"
            processes_to_use = 2 
    stem_cell_H1_HESC
    '''
    info_path_list = [("C:/Users/HP/Downloads/Siddharth/Data/hg38/GM12878/100kb/raw/chr1.mat.txt",(143000001,158000001),100000,10,"C:/Users/HP/Downloads/Siddharth/Data/hg38/GM12878/100kb/Chromosome_1","C:/Users/HP/Downloads/Siddharth/Thesis_Figures/Weighted_Network/GM12878/100kb/Chromosome_1"),
    ("C:/Users/HP/Downloads/Siddharth/Data/hg38/GM12878/10kb/raw/chr1.mat.txt",(143000001,158000001),10000,100,"C:/Users/HP/Downloads/Siddharth/Data/hg38/GM12878/10kb/Chromosome_1","C:/Users/HP/Downloads/Siddharth/Thesis_Figures/Weighted_Network/GM12878/10kb/Chromosome_1"),
    ("C:/Users/HP/Downloads/Siddharth/Data/hg38/stem_cell_H1_HESC/100kb/raw/chr1.mat.txt",(143000001,158000001),100000,10,"C:/Users/HP/Downloads/Siddharth/Data/hg38/stem_cell_H1_HESC/100kb/Chromosome_1","C:/Users/HP/Downloads/Siddharth/Thesis_Figures/Weighted_Network/stem_cell_H1_HESC/100kb/Chromosome_1"),
    ("C:/Users/HP/Downloads/Siddharth/Data/hg38/stem_cell_H1_HESC/10kb/raw/chr1.mat.txt",(143000001,158000001),10000,100,"C:/Users/HP/Downloads/Siddharth/Data/hg38/stem_cell_H1_HESC/10kb/Chromosome_1","C:/Users/HP/Downloads/Siddharth/Thesis_Figures/Weighted_Network/stem_cell_H1_HESC/10kb/Chromosome_1"),
    ("C:/Users/HP/Downloads/Siddharth/Data/hg38/GM12878/100kb/raw/chr2.mat.txt",(143000001,158000001),100000,10,"C:/Users/HP/Downloads/Siddharth/Data/hg38/GM12878/100kb/Chromosome_2","C:/Users/HP/Downloads/Siddharth/Thesis_Figures/Weighted_Network/GM12878/100kb/Chromosome_2"),
    ("C:/Users/HP/Downloads/Siddharth/Data/hg38/GM12878/10kb/raw/chr2.mat.txt",(143000001,158000001),10000,100,"C:/Users/HP/Downloads/Siddharth/Data/hg38/GM12878/10kb/Chromosome_2","C:/Users/HP/Downloads/Siddharth/Thesis_Figures/Weighted_Network/GM12878/10kb/Chromosome_2"),
    ("C:/Users/HP/Downloads/Siddharth/Data/hg38/stem_cell_H1_HESC/100kb/raw/chr2.mat.txt",(143000001,158000001),100000,10,"C:/Users/HP/Downloads/Siddharth/Data/hg38/stem_cell_H1_HESC/100kb/Chromosome_2","C:/Users/HP/Downloads/Siddharth/Thesis_Figures/Weighted_Network/stem_cell_H1_HESC/100kb/Chromosome_2"),
    ("C:/Users/HP/Downloads/Siddharth/Data/hg38/stem_cell_H1_HESC/10kb/raw/chr2.mat.txt",(143000001,158000001),10000,100,"C:/Users/HP/Downloads/Siddharth/Data/hg38/stem_cell_H1_HESC/10kb/Chromosome_2","C:/Users/HP/Downloads/Siddharth/Thesis_Figures/Weighted_Network/stem_cell_H1_HESC/10kb/Chromosome_2"),
    ("C:/Users/HP/Downloads/Siddharth/Data/hg38/GM12878/100kb/raw/chr3.mat.txt",(143000001,158000001),100000,10,"C:/Users/HP/Downloads/Siddharth/Data/hg38/GM12878/100kb/Chromosome_3","C:/Users/HP/Downloads/Siddharth/Thesis_Figures/Weighted_Network/GM12878/100kb/Chromosome_3"),
    ("C:/Users/HP/Downloads/Siddharth/Data/hg38/GM12878/10kb/raw/chr3.mat.txt",(143000001,158000001),10000,100,"C:/Users/HP/Downloads/Siddharth/Data/hg38/GM12878/10kb/Chromosome_3","C:/Users/HP/Downloads/Siddharth/Thesis_Figures/Weighted_Network/GM12878/10kb/Chromosome_3"),
    ("C:/Users/HP/Downloads/Siddharth/Data/hg38/stem_cell_H1_HESC/100kb/raw/chr3.mat.txt",(143000001,158000001),100000,10,"C:/Users/HP/Downloads/Siddharth/Data/hg38/stem_cell_H1_HESC/100kb/Chromosome_3","C:/Users/HP/Downloads/Siddharth/Thesis_Figures/Weighted_Network/stem_cell_H1_HESC/100kb/Chromosome_3"),
    ("C:/Users/HP/Downloads/Siddharth/Data/hg38/stem_cell_H1_HESC/10kb/raw/chr3.mat.txt",(143000001,158000001),10000,100,"C:/Users/HP/Downloads/Siddharth/Data/hg38/stem_cell_H1_HESC/10kb/Chromosome_3","C:/Users/HP/Downloads/Siddharth/Thesis_Figures/Weighted_Network/stem_cell_H1_HESC/10kb/Chromosome_3"),
    ("C:/Users/HP/Downloads/Siddharth/Data/hg38/GM12878/100kb/raw/chr7.mat.txt",(93140438,108140438),100000,10,"C:/Users/HP/Downloads/Siddharth/Data/hg38/GM12878/100kb/Chromosome_7","C:/Users/HP/Downloads/Siddharth/Thesis_Figures/Weighted_Network/GM12878/100kb/Chromosome_7"),
    ("C:/Users/HP/Downloads/Siddharth/Data/hg38/GM12878/10kb/raw/chr7.mat.txt",(93140438,108140438),10000,100,"C:/Users/HP/Downloads/Siddharth/Data/hg38/GM12878/10kb/Chromosome_7","C:/Users/HP/Downloads/Siddharth/Thesis_Figures/Weighted_Network/GM12878/10kb/Chromosome_7"),
    ("C:/Users/HP/Downloads/Siddharth/Data/hg38/stem_cell_H1_HESC/100kb/raw/chr7.mat.txt",(93140438,108140438),100000,10,"C:/Users/HP/Downloads/Siddharth/Data/hg38/stem_cell_H1_HESC/100kb/Chromosome_7","C:/Users/HP/Downloads/Siddharth/Thesis_Figures/Weighted_Network/stem_cell_H1_HESC/100kb/Chromosome_7"),
    ("C:/Users/HP/Downloads/Siddharth/Data/hg38/stem_cell_H1_HESC/10kb/raw/chr7.mat.txt",(93140438,108140438),10000,100,"C:/Users/HP/Downloads/Siddharth/Data/hg38/stem_cell_H1_HESC/10kb/Chromosome_7","C:/Users/HP/Downloads/Siddharth/Thesis_Figures/Weighted_Network/stem_cell_H1_HESC/10kb/Chromosome_7")]
    for info in info_path_list:
        HiC_data = matrix_reader(info[0])
        np.fill_diagonal(HiC_data,0) #Making diagonals of our matrix zero
        Regions_to_check = info[1]
        Resolution = info[2]
        node_no = info[3]
        data_root = info[4]
        image_root = info[5]
        processes_to_use = 2
        bins_to_check = (Regions_to_check[0]//Resolution,Regions_to_check[1]//Resolution)
        HiC_data_to_work = HiC_data[bins_to_check[0]:bins_to_check[1],bins_to_check[0]:bins_to_check[1]]
        Matrix_Forms_to_use_list = ["Raw","ICE","Z"]
        for Matrix_Form in Matrix_Forms_to_use_list:
            HiC_sub = multi_matrix_k_r_parallel_part_1(HiC_data_to_work,image_root,data_root,Matrix_Form,node_no,1)
            r_super_storer_1 = []
            k_super_storer_1 = []
            r_super_storer_2 = []
            k_super_storer_2 = []
            submatrix_sum_list = []
            max_ = len(HiC_sub)
            with Pool(processes=processes_to_use) as p, tqdm(total=max_) as pbar:
                    for result in p.starmap(k_r_calculator_parallel, zip((HiC_sub),repeat(node_no))):
                        pbar.update()
                        pbar.refresh()
                        r_super_storer_1.append(result[0])
                        k_super_storer_1.append(result[1])
                        r_super_storer_2.append(result[2])
                        k_super_storer_2.append(result[3])
                        submatrix_sum_list.append(result[4])
            matrix_writer("{}/{}/{}_k_values_kura_order".format(data_root,Matrix_Form,Matrix_Form),k_super_storer_1)
            matrix_writer("{}/{}/{}_r_values_kura_order".format(data_root,Matrix_Form,Matrix_Form),r_super_storer_1)
            matrix_writer("{}/{}/{}_k_values_uni_order".format(data_root,Matrix_Form,Matrix_Form),k_super_storer_2)
            matrix_writer("{}/{}/{}_r_values_uni_order".format(data_root,Matrix_Form,Matrix_Form),r_super_storer_2)
            list_writer("{}/{}/{}_submatrix_sum_list".format(data_root,Matrix_Form,Matrix_Form),submatrix_sum_list)

    '''
    ICE_normed = normalization.ICE_normalization(HiC_data_to_work)
    HiC_normed = submatrix_reader(ICE_normed,(node_no,node_no))
    r_super_storer_1 = []
    k_super_storer_1 = []
    r_super_storer_2 = []
    k_super_storer_2 = []
    submatrix_sum_list = []

    max_ = len(HiC_normed)
    with Pool(processes=2)as p, tqdm(total=max_) as pbar:
            for result in p.starmap(k_r_calculator_parallel, zip((HiC_normed),repeat(node_no))):
                pbar.update()
                pbar.refresh()
                r_super_storer_1.append(result[0])
                k_super_storer_1.append(result[1])
                r_super_storer_2.append(result[2])
                k_super_storer_2.append(result[3])
                submatrix_sum_list.append(result[4])
    matrix_writer("C:/Users/HP/Downloads/Siddharth/Data/hg38/GM12878/100kb/Chromosome_1/ICE_normalized/k_values_kura_order",k_super_storer_1)
    matrix_writer("C:/Users/HP/Downloads/Siddharth/Data/hg38/GM12878/100kb/Chromosome_1/ICE_normalized/r_values_kura_order",r_super_storer_1)
    matrix_writer("C:/Users/HP/Downloads/Siddharth/Data/hg38/GM12878/100kb/Chromosome_1/ICE_normalized/k_values_uni_order",k_super_storer_2)
    matrix_writer("C:/Users/HP/Downloads/Siddharth/Data/hg38/GM12878/100kb/Chromosome_1/ICE_normalized/r_values_uni_order",r_super_storer_2)
    list_writer("C:/Users/HP/Downloads/Siddharth/Data/hg38/GM12878/100kb/Chromosome_1/ICE_normalized/submatrix_sum_list",submatrix_sum_list)
    #print(results_list[0][0],"\n\n",len(results_list[0][0]),"\n\n",results_list[0][1],"\n\n",len(results_list[0][1]),results_list[0][2],"\n\n",len(results_list[0][2]),"\n\n",results_list[0][3],"\n\n",len(results_list[0][3]))
    time.sleep(60)
    z_normed = z_normalization(HiC_data)
    HiC_z_normed = submatrix_reader(z_normed,(node_no,node_no))
    r_super_storer_1 = []
    k_super_storer_1 = []
    r_super_storer_2 = []
    k_super_storer_2 = []
    submatrix_sum_list = []
    max_ = len(HiC_z_normed)
    with Pool(processes=2) as p, tqdm(total=max_) as pbar:
            for result in p.starmap(k_r_calculator_parallel, zip((HiC_z_normed),repeat(node_no))):
                pbar.update()
                pbar.refresh()
                r_super_storer_1.append(result[0])
                k_super_storer_1.append(result[1])
                r_super_storer_2.append(result[2])
                k_super_storer_2.append(result[3])
                submatrix_sum_list.append(result[4])
    matrix_writer("C:/Users/HP/Downloads/Siddharth/Data/hg38/GM12878/100kb/Chromosome_1/Z_normalized/k_values_kura_order",k_super_storer_1)
    matrix_writer("C:/Users/HP/Downloads/Siddharth/Data/hg38/GM12878/100kb/Chromosome_1/Z_normalized/r_values_kura_order",r_super_storer_1)
    matrix_writer("C:/Users/HP/Downloads/Siddharth/Data/hg38/GM12878/100kb/Chromosome_1/Z_normalized/k_values_uni_order",k_super_storer_2)
    matrix_writer("C:/Users/HP/Downloads/Siddharth/Data/hg38/GM12878/100kb/Chromosome_1/Z_normalized/r_values_uni_order",r_super_storer_2)
    list_writer("C:/Users/HP/Downloads/Siddharth/Data/hg38/GM12878/100kb/Chromosome_1/Z_normalized/submatrix_sum_list",submatrix_sum_list)
    '''