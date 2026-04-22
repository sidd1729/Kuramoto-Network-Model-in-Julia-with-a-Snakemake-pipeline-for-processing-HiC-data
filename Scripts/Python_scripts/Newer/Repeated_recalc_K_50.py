import numpy as np
from kuramoto import Kuramoto
from tqdm import tqdm
from iced import normalization
import multiprocessing
from multiprocessing import Pool
from itertools import repeat
import csv
import pickle
import pandas as pd
import pyarrow.feather as feather
import logging
import sys
#Efficiency improvement:
def derivative(self, angles_vec, t, adj_mat, coupling):
        '''
        Compute derivative of all nodes for current state, defined as

        dx_i    natfreq_i + k  sum_j ( Aij* sin (angle_j - angle_i) )
        ---- =             ---
         dt                M_i

        t: for compatibility with scipy.odeint
        '''
        assert len(angles_vec) == len(self.natfreqs) == len(adj_mat), \
            'Input dimensions do not match, check lengths'
        interactions = adj_mat * np.sin(angles_vec[:,None] - angles_vec)  # Aij * sin(j-i)

        dxdt = self.natfreqs + coupling * interactions.sum(axis=0)  # sum over incoming interactions
        return dxdt
Kuramoto.derivative = derivative

#Pickling for ease of access, and don't need to worry about how to write to storage:
def general_pickler(filename,tbp_object=None,read = False,reader_fn = None):
    if read == True:
        tbp_object = reader_fn(filename)      
    with open(filename + ".pickle",'wb') as f:
        pickle.dump(tbp_object,f)

def submatrix_reader(matrix,submatrix_shape,stride = 1,offset = (0,0),pattern = "diagonal"):
    #Logic should be working only for symmetric submatrices, not properly generalized to non square matrices
    submatrix_no = min(matrix.shape) - max(submatrix_shape) + 1 - max(offset)
    submatrices=[]
    for i in range(0,submatrix_no,stride):
        submatrix = matrix[i+offset[0]:i+submatrix_shape[0]+offset[0],i+offset[1]:i+submatrix_shape[1] + offset[1]]
        submatrices.append(submatrix)
    return submatrices

def matrix_reader(filename,delete = False):
    if filename[-3:] == "npy":
        return np.load(filename)
    elif filename[-7:] == "feather":
        np_mat = feather.read_feather(filename).to_numpy()
        return np_mat
    else:
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
    
def z_normalization(matrix,iteration_no = 10):
    #Numpy arrays are mutable so we first make a copy to preserve the original array and then use mutability to our advantage
    matrix = matrix.copy()
    for _ in range(0,iteration_no):
        for i in range(0,matrix.shape[0]):
            matrix[i,:] = (matrix[i,:]- np.mean(matrix[i,:]))/(np.std(matrix[i,:])+1e-48)
        for j in range(0,matrix.shape[1]):
            matrix[:,j] = (matrix[:,j]- np.mean(matrix[:,j]))/(np.std(matrix[:,j])+1e-48)
    return matrix

def tuple_list_writer(filename,variable_name):
    file = open('{}'.format(filename), 'w')
    for t in variable_name:
      file.write(' '.join(str(s) for s in t) + '\n')
    file.close()

def universal_order_parameter_modif_2(phase_time_series_matrix,adj_mat):
        '''
        Compute universal order parameter R_t - mean length of resultant vector
        '''
        #suma = 0
        
        angles_vec = phase_time_series_matrix
        #Broadcasting: Mind blowing!?
        angles_i_j_sum = (np.cos(angles_vec[:,np.newaxis,:] - angles_vec)).sum(-1)
        #print(angles_i_j_sum.shape)
        #print(angles_i.shape,angles_j.shape)
        angles_i_j = np.divide(angles_i_j_sum,angles_vec.shape[1])
        suma = np.sum(adj_mat*(angles_i_j))
    #print("{} is overall sum, \n {} is added sum, \n {} is adjacency mat weight, \n {} is difference between angle vectors".format(suma,adj_mat[i,j]*sum(np.cos(i_vec[-1000:]-j_vec[-1000:]))/len(i_vec[-1000:]),adj_mat[i,j],i_vec[-1000:]-j_vec[-1000:]))
    #else:
    #suma_term = adj_mat[i,j]*np.sum(np.cos(i_vec-j_vec))/len(i_vec)
    #return suma_term
        #suma = 0                
        denom = np.sum(adj_mat)
        ##suma = functools.reduce(lambda x, y: x+y, suma_appl_list)
        #suma = sum(suma_appl_list)
        #if abs(suma/denom) > 1:
            #return (abs(suma/denom),phase_time_series_matrix,denom)
        #print(suma,denom)
        return (abs(suma / denom),0,denom) #Change back after trials

def K_50_recalc_init(HiC_data_file,Regions_to_check,Resolution,node_no,data_root,Matrix_Form,Order_parameter_form,to_write = False):
    HiC_data = matrix_reader(HiC_data_file)
    np.fill_diagonal(HiC_data,0) #Making diagonals of our matrix zero
    bins_to_check = (Regions_to_check[0]//Resolution,Regions_to_check[1]//Resolution)
    #print(bins_to_check)
    HiC_data_to_work = HiC_data[bins_to_check[0]:bins_to_check[1],bins_to_check[0]:bins_to_check[1]]
    if Matrix_Form == "Raw":
        print("Hit Raw")
        HiC_data_to_work = HiC_data_to_work
    elif Matrix_Form == "ICE":
        print("Hit ICE")
        HiC_data_to_work = normalization.ICE_normalization(HiC_data_to_work)
    elif Matrix_Form == "Z":
        print("Hit Z")
        HiC_data_to_work = z_normalization(HiC_data_to_work)
    HiC_sub = submatrix_reader(HiC_data_to_work,(node_no,node_no))
    #yan_val_storer_1 = tuple_list_reader("{}/{}/{}/{}_K_50".format(data_root,Matrix_Form,Order_parameter_form,Order_parameter_form)) #### To be changed back when code works properly
    #empty_HiC_list = list_reader("{}/{}/{}/{}_zero_index_list".format(data_root,Matrix_Form,Order_parameter_form,Order_parameter_form))
    #r_1_list_1 = tuple_list_reader("{}/{}/{}/{}_r_1".format(data_root,Matrix_Form,Order_parameter_form,Order_parameter_form))
    #non_positive_roots_1 = tuple_list_reader("{}/{}/{}/{}_non_pos_roots".format(data_root,Matrix_Form,Order_parameter_form,Order_parameter_form))
    #K_50_indices_list = list(map(lambda x : x[1],yan_val_storer_1)) #### To be changed back when code works properly
    K_50_indices_list = []
    Possible_indices_list = list(range(0,len(HiC_sub)))
    print(Possible_indices_list)
    indices_to_act_list = list(set(Possible_indices_list) - set(K_50_indices_list))
    return (indices_to_act_list,HiC_sub,node_no)

def k_r_calculator_parallel(HiC_submatrix_entry,no_nodes,coupling_range = (0,0.6),no_coupling=200,adaptive_coupling=False,adaptive_region = (0,0.6),adaptive_cut = 100):
        if adaptive_coupling == True:
            coupling_vals_part_1 = np.linspace(adaptive_region[0],adaptive_region[1],adaptive_cut)
            coupling_vals_part_2 = np.linspace(adaptive_region[1],coupling_range[1],no_coupling-adaptive_cut)
            coupling_vals = np.unique(np.hstack((coupling_vals_part_1,coupling_vals_part_2)))
            #print("Coupled")
        else:
            coupling_vals = np.linspace(coupling_range[0],coupling_range[1],no_coupling)
        n_nodes=no_nodes
        coupling = coupling_vals
        #n_nodes=10
        #graph = HiC_50[val]
        graph = HiC_submatrix_entry
        adjmat=graph
        runs = []
        
        for coupling in tqdm(coupling_vals):
            model = Kuramoto(coupling=coupling, dt=0.1, T=500, n_nodes=n_nodes) 
            model.natfreqs = np.random.normal(1, 0.1, size=n_nodes)  # reset natural frequencies
            act_mat = model.run(adj_mat=adjmat)
            runs.append(act_mat)
        runs_array = np.array(runs)
        r_k_storer_1 = []
        r_k_storer_2 = []
        for i, coupling in tqdm(enumerate(coupling_vals)):
            r_mean_kura = np.mean([model.phase_coherence(vec)
                          for vec in runs_array[i, :, -1000:].T]) # mean over last 1000 steps
            r_mean_uni_modif_2 = universal_order_parameter_modif_2(runs_array[i],adjmat)[0]
            r_mean_uni = r_mean_uni_modif_2
            
            r_k_storer_1.append((coupling,r_mean_kura))
            r_k_storer_2.append((coupling,r_mean_uni))
        r_storer_1 = list(map(lambda x: x[1],r_k_storer_1))
        k_storer_1 = list(map(lambda x: x[0],r_k_storer_1))

        r_storer_2 = list(map(lambda x: x[1],r_k_storer_2))
        k_storer_2 = list(map(lambda x: x[0],r_k_storer_2))
        return (r_storer_1,k_storer_1,r_storer_2,k_storer_2,np.sum(HiC_submatrix_entry))



def K_50_fit_checker(indices_to_act_list,HiC_sub,no_nodes,repetition_rate,storage_path_root = None,r_val_filename = None ,k_val_filename = None):
    changed_root_list =[]
    init_interval = (0.0,0.6)
    for index in indices_to_act_list:
        test_entry = HiC_sub[indices_to_act_list[index]] #10 was the index previously used for the other tests, where hyperbolic tan fit very well, #5th bin for 100kb and 50th for 10kb
        no_coupling = 500
        #print(init_partitions)
        r_super_storer = []
        k_super_storer = []
        for _ in range(repetition_rate):
            r_storer_1,k_storer_1,r_storer_2,k_storer_2,submat_sum = k_r_calculator_parallel(test_entry,no_nodes,coupling_range = init_interval,no_coupling=no_coupling,adaptive_coupling= True,adaptive_region = (0,0.025),adaptive_cut=500)
            r_super_storer.append(r_storer_1)
            k_super_storer.append(k_storer_1)
        r_super_array = np.array(r_super_storer)
        k_super_array = np.array(k_super_storer)
        mean_stddev_list = []
        complementary_mean_stddev_list = []
        for j in range(0,k_super_array.shape[1]):   
            entry_slice = r_super_array[:,j]
            entry_mean = np.mean(entry_slice)
            entry_stdev = np.std(entry_slice)
            mean_stddev_list.append((entry_mean,entry_stdev))
            complementary_slice = k_super_array[:,j]
            complementary_mean = np.mean(complementary_slice)
            complementary_stdev = np.std(complementary_slice)
            complementary_mean_stddev_list.append((complementary_mean,complementary_stdev))
        #tuple_list_writer("./GM12878/ICE_5_5/10k_coupling_trial_500_Chr1_100kb__GM12878_ICE_5_5_first_15_{}".format(index),mean_stddev_list)
        #tuple_list_writer("./GM12878/ICE_5_5/10k_coupling_k_vals__500_Chr1_100kb_GM12878_ICE_5_5_first_15_{}".format(index),complementary_mean_stddev_list)
    return changed_root_list

def K_50_fit_checker_parallel(HiC_entry,no_nodes,repetition_rate = 100,init_interval = (0,0.6),no_coupling = 500,adaptive_coupling = True,adaptive_region = (0,0.025),adaptive_cut = 500):
     r_super_storer_1 = []
     k_super_storer_1 = []
     r_super_storer_2 = []
     k_super_storer_2 = []
     for _ in range(repetition_rate):
            r_storer_1,k_storer_1,r_storer_2,k_storer_2,submat_sum = k_r_calculator_parallel(HiC_entry,no_nodes,coupling_range = init_interval,no_coupling=no_coupling,adaptive_coupling= adaptive_coupling,adaptive_region = adaptive_region,adaptive_cut=adaptive_cut)
            r_super_storer_1.append(r_storer_1)
            k_super_storer_1.append(k_storer_1)
            r_super_storer_2.append(r_storer_2)
            k_super_storer_2.append(k_storer_2)
     r_super_array_1 = np.array(r_super_storer_1)
     k_super_array_1 = np.array(k_super_storer_1)
     r_super_array_2 = np.array(r_super_storer_2)
     k_super_array_2 = np.array(k_super_storer_2)
     mean_stddev_list = []
     complementary_mean_stddev_list = []
     mean_stddev_list_2 = []
     complementary_mean_stddev_list_2 = []
     for j in range(0,k_super_array_1.shape[1]):   
            entry_slice = r_super_array_1[:,j]
            entry_mean = np.mean(entry_slice)
            entry_stdev = np.std(entry_slice)
            mean_stddev_list.append((entry_mean,entry_stdev))
            complementary_slice = k_super_array_1[:,j]
            complementary_mean = np.mean(complementary_slice)
            complementary_stdev = np.std(complementary_slice)
            complementary_mean_stddev_list.append((complementary_mean,complementary_stdev))
     
            entry_slice_2 = r_super_array_2[:,j]
            entry_mean_2 = np.mean(entry_slice_2)
            entry_stdev_2 = np.std(entry_slice_2)
            mean_stddev_list_2.append((entry_mean_2,entry_stdev_2))
            complementary_slice_2 = k_super_array_2[:,j]
            complementary_mean_2 = np.mean(complementary_slice_2)
            complementary_stdev_2 = np.std(complementary_slice_2)
            complementary_mean_stddev_list_2.append((complementary_mean_2,complementary_stdev_2))
     
     return (mean_stddev_list,complementary_mean_stddev_list,mean_stddev_list_2,complementary_mean_stddev_list_2)#,HiC_entry)
     
     
     

if __name__ == '__main__':
    #multiprocessing.set_start_method('spawn')
    supreme_k_storer_1 = []
    supreme_k_storer_2 = []
    supreme_r_storer_1 = []
    supreme_r_storer_2 = []                    #/home/ksslab/Siddharth/Program_Directory/Data/hg38/GM12878/100kb/raw/chr1.mat.txt       #(143000001,158000001)                                                                              #(0,248956425)
    test,test_mat, no_nodes = K_50_recalc_init(sys.argv[1],(15000000,30000000),40000,5,"/home/ksslab/Siddharth/Program_Directory/Data/hg38/GM12878/100kb/Chromosome_1","ICE","Kuramoto")
    print(test)
    processes_to_use  = 24
    max_ = len(test_mat)
    #max_ = 2
    #counter_var = 0
    with Pool(processes=processes_to_use) as p, tqdm(total=max_) as pbar:
                    for result in p.starmap(K_50_fit_checker_parallel,zip(test_mat,repeat(no_nodes),repeat(100))):
                        pbar.update()
                        pbar.refresh()
                        #counter_var+=1
                        #print("Success !,{}".format(counter_var))
                        supreme_r_storer_1.append(result[0])
                        supreme_k_storer_1.append(result[1])
                        supreme_r_storer_2.append(result[2])
                        supreme_k_storer_2.append(result[3])
                        #print(np.where(test_mat == result[4]))
    #feather.write_feather(pd.DataFrame(supreme_r_storer_1),"/home/ksslab/Siddharth/Program_Directory/Data/Simulation_Results/Linux/GM12878/143_158_Reproducibility/1000_Repeat/Chr1_Kura_r_1.feather",)
    #feather.write_feather(pd.DataFrame(supreme_k_storer_1),"/home/ksslab/Siddharth/Program_Directory/Data/Simulation_Results/Linux/GM12878/143_158_Reproducibility/1000_Repeat/Chr1_Kura_k_1.feather")
    #feather.write_feather(pd.DataFrame(supreme_r_storer_2),"/home/ksslab/Siddharth/Program_Directory/Data/Simulation_Results/Linux/GM12878/143_158_Reproducibility/1000_Repeat/Chr1_Uni_r_1.feather")
    #feather.write_feather(pd.DataFrame(supreme_k_storer_2),"/home/ksslab/Siddharth/Program_Directory/Data/Simulation_Results/Linux/GM12878/143_158_Reproducibility/1000_Repeat/Chr1_Uni_k_1.feather")
    feather.write_feather(pd.DataFrame(supreme_r_storer_1),sys.argv[2])
    feather.write_feather(pd.DataFrame(supreme_k_storer_1),sys.argv[3])
    feather.write_feather(pd.DataFrame(supreme_r_storer_2),sys.argv[4])
    feather.write_feather(pd.DataFrame(supreme_k_storer_2),sys.argv[5])