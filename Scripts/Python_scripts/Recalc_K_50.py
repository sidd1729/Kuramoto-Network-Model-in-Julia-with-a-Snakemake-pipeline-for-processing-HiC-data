import matplotlib.pyplot as plt
import numpy as np
from iced import normalization
from tqdm import tqdm
import time
import csv
from kuramoto import Kuramoto
import os
import scipy
#import profile
#We check the 4 lists and if K_50 is absent, we divide our interval into n subintervals, for every subinterval, we calculate minimum and maximum, we can then select those
#subintervals, where r = 0.5 for some K value is likely, can we necessarily make a spline for this interval?

#Thus for non_positive and r_1 values, we begin this partition or actually non K_50 values, how to decide convergence, somehow need recurrent evaluation of partition at one level and repeat recurrent evaluation till we reach convergence of K_0.5

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

def submatrix_reader(matrix,submatrix_shape,stride = 1,offset = (0,0),pattern = "diagonal"):
    #Logic should be working only for symmetric submatrices, not properly generalized to non square matrices
    submatrix_no = min(matrix.shape) - max(submatrix_shape) + 1 - max(offset)
    submatrices=[]
    for i in range(0,submatrix_no,stride):
        submatrix = matrix[i+offset[0]:i+submatrix_shape[0]+offset[0],i+offset[1]:i+submatrix_shape[1] + offset[1]]
        submatrices.append(submatrix)
    return submatrices

def z_normalization(matrix,iteration_no = 10):
    #Numpy arrays are mutable so we first make a copy to preserve the original array and then use mutability to our advantage
    matrix = matrix.copy()
    for _ in range(0,iteration_no):
        for i in range(0,matrix.shape[0]):
            matrix[i,:] = (matrix[i,:]- np.mean(matrix[i,:]))/(np.std(matrix[i,:])+1e-48)
        for j in range(0,matrix.shape[1]):
            matrix[:,j] = (matrix[:,j]- np.mean(matrix[:,j]))/(np.std(matrix[:,j])+1e-48)
    return matrix

def list_reader(filename):
    file = open('{}'.format(filename), 'r')
    content = file.readlines()
    new_item_list =[]
    #print(len(content[1]))
    file.close()
    tsv_reader = csv.reader(content, delimiter='\t')
    for item_list in tsv_reader:
        #del(item_list[len(item_list)-1]) #Only to be used when there is extra blank character at the end
        for item in item_list:
            new_item_list.append(float(item))
    return new_item_list

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

def tuple_list_writer(filename,variable_name):
    file = open('{}'.format(filename), 'w')
    for t in variable_name:
      file.write(' '.join(str(s) for s in t) + '\n')
    file.close()

def tuple_list_reader(filename):
    file = open('{}'.format(filename), 'r')
    content = file.read().splitlines()
    tuple_list =[]
    file.close()
    for row in content:
        tuple_list.append(tuple(float(i) for i in row.split(' ')))
    return tuple_list

def Region_partition(region_start,region_stop,partition_no):
    new_partitions_points_list = []
    new_partitions_list = [0]*partition_no
    new_partitions_list_alt  = [0]*partition_no
    partition_points_alt = np.linspace(region_start,region_stop,partition_no+2)
    for i in range(1,partition_no):
        partition_points = (region_start+region_stop)* i/partition_no
        new_partitions_points_list.append(partition_points)
    for i in range(partition_no):
        if i == 0:
            new_partitions_list[i] = (region_start,new_partitions_points_list[i])
            new_partitions_list_alt[i] = (region_start,partition_points_alt[i+1])
        elif i == partition_no-1:
            new_partitions_list[i] = (new_partitions_points_list[i-1],region_stop)
            new_partitions_list_alt[i] = (partition_points_alt[i],region_stop)
        else:
            new_partitions_list[i] = (new_partitions_points_list[i-1] , new_partitions_points_list[i])
            new_partitions_list_alt[i] = (partition_points_alt[i],partition_points_alt[i+1])
    return (new_partitions_list,new_partitions_points_list,partition_points_alt,new_partitions_list_alt)

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

def k_r_calculator_parallel(HiC_submatrix_entry,no_nodes,coupling_range = (0,0.6),no_coupling=200,adaptive_coupling=False,adaptive_region = (0,0.6),adaptive_cut = 100):
        #k_50_storer_1 = []
        #k_50_storer_2 =[]
        if adaptive_coupling == True:
            coupling_vals_part_1 = np.linspace(adaptive_region[0],adaptive_region[1],adaptive_cut)
            coupling_vals_part_2 = np.linspace(adaptive_region[1],coupling_range[1],no_coupling-adaptive_cut)
            coupling_vals = np.unique(np.hstack((coupling_vals_part_1,coupling_vals_part_2)))
            #print("Coupled")
        else:
            coupling_vals = np.linspace(coupling_range[0],coupling_range[1],no_coupling)
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
        #plt.figure()
        #plt.hist(model.natfreqs)
        #plt.xlabel('natural frequency')
        #plt.ylabel('count')
        #start_time = time.perf_counter()
        # Plot all time series for all coupling values (color coded)
        runs_array = np.array(runs)
        #end_time = time.perf_counter()
        #elapsed_time = end_time - start_time
        #print(f"Execution time: {elapsed_time:.6f} seconds")
        #plt.show()
        #plt.close()
        r_k_storer_1 = []
        r_k_storer_2 = []
        #plt.figure()
        for i, coupling in tqdm(enumerate(coupling_vals)):
            r_mean_kura = np.mean([model.phase_coherence(vec)
                          for vec in runs_array[i, :, -1000:].T]) # mean over last 1000 steps
            r_mean_uni = universal_order_parameter(runs_array[i],adjmat)
            r_k_storer_1.append((coupling,r_mean_kura))
            r_k_storer_2.append((coupling,r_mean_uni))
            #plt.scatter(coupling, r_mean_kura, c='steelblue', s=20, alpha=0.7,label = "Kuramoto")
            #plt.scatter(coupling, r_mean_uni, c='red', s=20, alpha=0.7,label = "Universal")

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
        #plt.close()
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

def k_50_root_routine(k_vals_list,r_vals_list):
        list_writer("./r_trial",r_vals_list)
        list_writer("./k_trial",k_vals_list)
        #print("The k and r values for the problem are \n {} and  \n\n   {} \n\n".format(k_vals_list,r_vals_list))
        #Not done yet remember to transform the values back/implement an inverse transform and pruning
        try:
            spl = scipy.interpolate.make_smoothing_spline(k_vals_list,r_vals_list)
        except:
            k_map_min_max_transformed = []
            for i in range(0,len(k_vals_list)):
                k_map_min_max_transformed.append((k_vals_list[i]-min(k_vals_list))/(max(k_vals_list) - min(k_vals_list)))
            inverse_transformed = []
            for i in range(0,len(k_vals_list)):
                inverse_transformed.append(k_map_min_max_transformed[i]*(max(k_vals_list) - min(k_vals_list)) + min(k_vals_list))
            #print("The k and r values for the ill-formed problem are \n {} and {} \n\n".format(k_vals_list[0:10],r_vals_list[0:10]))
            #print("The min-max transformed value is {}".format(k_map_min_max_transformed))
            #print("Our inverse transform is {}".format(inverse_transformed[0:10]))
            #try:
            spl = scipy.interpolate.make_smoothing_spline(k_map_min_max_transformed,r_vals_list)
            #except:
        #try:        
        #spl_0_5 = scipy.interpolate.make_smoothing_spline(k_vals_list,np.array(r_vals_list)-0.5)
        #xvals = np.linspace(0, 0.6, 200)
        try:
            spl_0_5 = scipy.interpolate.make_smoothing_spline(k_vals_list,np.array(r_vals_list)-0.5)
        except:
            k_map_min_max_transformed =[]
            for i in range(0,len(k_vals_list)):
                k_map_min_max_transformed.append((k_vals_list[i]-min(k_vals_list))/(max(k_vals_list) - min(k_vals_list)))
            inverse_transformed = []
            for i in range(0,len(k_vals_list)):
                inverse_transformed.append(k_map_min_max_transformed[i]*(max(k_vals_list) - min(k_vals_list)) + min(k_vals_list))
            #print("The k and r values for the ill-formed problem are \n {} and {} \n\n".format(k_vals_list,r_vals_list))
            #print("The min-max transformed value is {}".format(k_map_min_max_transformed))
           #print("Our inverse transform is {}".format(inverse_transformed[0:10]))
            spl_0_5 = scipy.interpolate.make_smoothing_spline(k_map_min_max_transformed,np.array(r_vals_list)-0.5)
        ppoly = scipy.interpolate.PPoly.from_spline(spl_0_5)
        #print(ppoly.x)
        #print(ppoly.c.T)
        #for i in range(ppoly.x):
           #poly_string = ''
           #for j in range(len(ppoly.c[:,i])):
                
        #evaluated_0_5_vals = spl(ppoly.roots())
        #xvals = np.linspace(0, 0.6, 200)
        if len(ppoly.roots(extrapolate = False)) >=3:
            filtered_root = [j for j in ppoly.roots(extrapolate = False) if j > 0 and j < max(ppoly.roots())]
        if len(ppoly.roots(extrapolate = False)) ==  2:
            filtered_root = [j for j in ppoly.roots(extrapolate = False ) if j > 0]
        elif len(ppoly.roots(extrapolate = False)) == 1:
            filtered_root = ppoly.roots(extrapolate = False)
        elif len(ppoly.roots(extrapolate = False)) == 0:
            filtered_root = []
        if len(filtered_root)>0:
            filtered_root = filtered_root[0]
        elif len(filtered_root) == 0:
            filtered_root = 2000
        else:
            filtered_root = ppoly.roots()[0]
        if 'k_map_min_max_transformed' in locals():
            return ((filtered_root*(max(k_vals_list) - min(k_vals_list)) + min(k_vals_list)),spl(filtered_root))
        return (filtered_root,spl(filtered_root))
def K_50_recalc_init(HiC_data_file,Regions_to_check,Resolution,node_no,data_root,Matrix_Form,Order_parameter_form,to_write = False):
    HiC_data = matrix_reader(HiC_data_file)
    np.fill_diagonal(HiC_data,0) #Making diagonals of our matrix zero
    bins_to_check = (Regions_to_check[0]//Resolution,Regions_to_check[1]//Resolution)
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

def K_50_recalc_manager(indices_to_act_list,HiC_sub,no_nodes,partition_rate):
    init_interval = (0.0,0.6)
    #for index in indices_to_act_list:
    test_entry = HiC_sub[indices_to_act_list[7]]
    init_partitions = (Region_partition(init_interval[0],init_interval[1],partition_rate))[3]
    #print(init_partitions)
    changed_root_list = []
    for partition in init_partitions:
        r_storer_1,k_storer_1,r_storer_2,k_storer_2,submat_sum = k_r_calculator_parallel(test_entry,no_nodes,coupling_range = partition,no_coupling=200)
        #print(set(r_storer_1))
        if (len(set(r_storer_1)) == 1):
            print("Always 1")
            continue
        #print(r_storer_1,k_storer_1)
        Nan_status_list = map(np.isnan,r_storer_1)
        ind = [i for i, val in enumerate(Nan_status_list) if val == True]
        #print(ind)
        for indice in ind:
            #print(r_storer_1[indice])
            if r_storer_1[indice-1] == r_storer_1[indice+1]:
                r_storer_1[indice] = r_storer_1[indice+1]
        filtered_root,eval_filtered_root = k_50_root_routine(k_storer_1,r_storer_1)
        if eval_filtered_root == 0.5 and filtered_root < 1 and filtered_root > 0:
            changed_entry = (filtered_root,eval_filtered_root)
            changed_root_list.append(changed_entry)
            break
    iteration_no = 0
    if changed_root_list == []:
        iteration_counter = 0 
        while iteration_counter < iteration_no and changed_root_list == []:
             all_partitions = []
             #if iteration_counter >=1:
                #init_partitions = subpartitions 
             print("Iteration no is {}".format(iteration_counter))#,init_partitions)
             for partition in init_partitions:
               subpartitions = (Region_partition(partition[0],partition[1],partition_rate))[3]
               all_partitions.extend(subpartitions)
               print("{} and {} are lengths of current subpartitions list and all subparitions list".format(len(subpartitions),len(all_partitions)))
               for subpartition in subpartitions:
                   print("{} is subpartition we are currently looking at ".format(subpartition))
                   r_storer_1,k_storer_1,r_storer_2,k_storer_2,submat_sum = k_r_calculator_parallel(test_entry,no_nodes,coupling_range = subpartition,no_coupling=200)
                   #if not min(r_storer_1) < 0.5 < max(r_storer_1):
                       
                   #print(set(r_storer_1))
                   if (len(set(r_storer_1)) == 1):
                        print("Always 1")
                        continue
                   #print(r_storer_1,k_storer_1)
                   #print(k_storer_1)
                   Nan_status_list = map(np.isnan,r_storer_1)
                   ind = [i for i, val in enumerate(Nan_status_list) if val == True]
                   #print(ind)
                   for indice in ind:
                       #print(r_storer_1[indice])
                       if r_storer_1[indice-1] == r_storer_1[indice+1]:
                           r_storer_1[indice] = r_storer_1[indice+1]
                   filtered_root,eval_filtered_root = k_50_root_routine(k_storer_1,r_storer_1)
                   if eval_filtered_root == 0.5 and filtered_root < 1 and filtered_root > 0:
                        changed_entry = (filtered_root,eval_filtered_root)
                        changed_root_list.append(changed_entry)       
               iteration_counter+=1
               init_partitions = all_partitions 
    return changed_root_list

def K_50_recalc_manager_sampling_vers(indices_to_act_list,HiC_sub,no_nodes,sampling_rate):
    init_interval = (0.0,0.6)
    #for index in indices_to_act_list:
    test_entry = HiC_sub[indices_to_act_list[5]] #10 was the index previously used for the other tests, where hyperbolic tan fit very well
    no_coupling = 100000
    #print(init_partitions)
    changed_root_list = []
    iteration_no = 1
    if changed_root_list == []:
        iteration_counter = 0 
        while iteration_counter < iteration_no and changed_root_list == []:
            no_coupling *= sampling_rate 
            r_storer_1,k_storer_1,r_storer_2,k_storer_2,submat_sum = k_r_calculator_parallel(test_entry,no_nodes,coupling_range = init_interval,no_coupling=no_coupling,adaptive_coupling= False,adaptive_region = (0,0.35),adaptive_cut=8000)
            if (len(set(r_storer_1)) == 1):
                print("Always 1")
                continue
            #print(len(r_storer_1),len(k_storer_1))
            #print(k_storer_1)
            Nan_status_list = map(np.isnan,r_storer_1)
            ind = [i for i, val in enumerate(Nan_status_list) if val == True]
            #print(ind)
            for indice in ind:
                #print(r_storer_1[indice])
                if r_storer_1[indice-1] == r_storer_1[indice+1]:
                    r_storer_1[indice] = r_storer_1[indice+1]
            filtered_root,eval_filtered_root = k_50_root_routine(k_storer_1,r_storer_1)
            if eval_filtered_root == 0.5 and filtered_root < 1 and filtered_root > 0:
                    changed_entry = (filtered_root,eval_filtered_root)
                    changed_root_list.append(changed_entry)       
            iteration_counter+=1
    return changed_root_list

def K_50_fit_checker(indices_to_act_list,HiC_sub,no_nodes,repetition_rate):
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
        tuple_list_writer("./GM12878/ICE_5_5/10k_coupling_trial_500_Chr1_100kb__GM12878_ICE_5_5_first_15_{}".format(index),mean_stddev_list)
        tuple_list_writer("./GM12878/ICE_5_5/10k_coupling_k_vals__500_Chr1_100kb_GM12878_ICE_5_5_first_15_{}".format(index),complementary_mean_stddev_list)
    return changed_root_list

test,test_mat, no_nodes = K_50_recalc_init("C:/Users/HP/Downloads/Siddharth/Data/hg38/GM12878/100kb/raw/chr1.mat.txt",(143000001,158000001),100000,5,"C:/Users/HP/Downloads/Siddharth/Data/hg38/GM12878/100kb/Chromosome_1","ICE","Kuramoto")
#temp_indices_to_act = [1,2,3,4,5,6,7,8,9,10,11]
#temp_indices_to_act_tuples = tuple_list_reader("C:/Users/HP/Downloads/Siddharth/Data/hg38/GM12878/100kb/Chromosome_1/ICE/Universal/Universal_K_50")
#temp_indices_to_act =  list(map(int,map(lambda x: x[1],temp_indices_to_act_tuples)))
#print(temp_indices_to_act,len(temp_indices_to_act))
print(test)
changed_root_list = K_50_fit_checker(test,test_mat,no_nodes,100)     #Fit checking
#changed_root_list = K_50_recalc_manager_sampling_vers(test,test_mat,no_nodes,1) #Actual Line
#changed_root_list = K_50_recalc_manager(test,test_mat,no_nodes,50)
#print(changed_root_list)
#partitions,partition_points,partition_points_alt,new_partitions_list_alt = Region_partition(0.12,0.24,5)
#print("{}\n{}\n{}\n{}".format(partitions,partition_points,partition_points_alt,new_partitions_list_alt))

'''
    if Order_parameter_form == "Kuramoto":
        k_super_storer_1 =  matrix_reader("{}/{}/{}_k_values_kura_order".format(data_root,Matrix_Form,Matrix_Form))
        r_super_storer_1 =  matrix_reader("{}/{}/{}_r_values_kura_order".format(data_root,Matrix_Form,Matrix_Form))
    elif Order_parameter_form == "Universal":
        k_super_storer_1 =  matrix_reader("{}/{}/{}_k_values_uni_order".format(data_root,Matrix_Form,Matrix_Form))
        r_super_storer_1 =  matrix_reader("{}/{}/{}_r_values_uni_order".format(data_root,Matrix_Form,Matrix_Form))
    yan_val_storer_1 = []
    empty_HiC_list =[]
    r_1_list_1 = []
    non_positive_roots_1 = []
    for i in tqdm(range(0,len(k_super_storer_1))):
        if np.sum(HiC_sub[i]) == 0:
            empty_HiC_list.append(i)
            continue
        x = k_super_storer_1[i]
        y = r_super_storer_1[i]
        spl = scipy.interpolate.make_smoothing_spline(x,y)
        xvals = np.linspace(0, 0.6, 200)
        spl_0_5 = scipy.interpolate.make_smoothing_spline(x,np.array(y)-0.5)
        ppoly = scipy.interpolate.PPoly.from_spline(spl_0_5)
        evaluated_0_5_vals = spl(ppoly.roots())
        xvals = np.linspace(0, 0.6, 200)
        if len(ppoly.roots()) >=3:
            filtered_root = [j for j in ppoly.roots() if j > 0 and j < max(ppoly.roots())]
        if len(ppoly.roots()) == 2:
            filtered_root = [j for j in ppoly.roots() if j > 0]
        elif len(ppoly.roots()) == 1:
            filtered_root = ppoly.roots()
        elif len(ppoly.roots()) == 0:
            filtered_root = []
        if len(filtered_root)>0:
            filtered_root = filtered_root[0]
        elif len(filtered_root) == 0:
            filtered_root = 2000
        else:
            filtered_root = ppoly.roots()[0]
        if filtered_root > 2:
            r_1_list_1.append((filtered_root,i))
        elif filtered_root < 0 :
            non_positive_roots_1.append((filtered_root,i))
        else:
            yan_val_storer_1.append((filtered_root,i))
    Data_Dir_exist = os.path.exists("{}/{}/{}".format(data_root,Matrix_Form,Order_parameter_form)) 
    #print("Prior to dir", Data_Dir_exist)                                          
    if not Data_Dir_exist:
         os.makedirs("{}/{}/{}".format(data_root,Matrix_Form,Order_parameter_form))
         print("Made")
    #print("Reached the other side")
    if to_write == True:
        tuple_list_writer("{}/{}/{}/{}_K_50".format(data_root,Matrix_Form,Order_parameter_form,Order_parameter_form),yan_val_storer_1)
        list_writer("{}/{}/{}/{}_zero_index_list".format(data_root,Matrix_Form,Order_parameter_form,Order_parameter_form),empty_HiC_list)
        tuple_list_writer("{}/{}/{}/{}_r_1".format(data_root,Matrix_Form,Order_parameter_form,Order_parameter_form),r_1_list_1)
        tuple_list_writer("{}/{}/{}/{}_non_pos_roots".format(data_root,Matrix_Form,Order_parameter_form,Order_parameter_form),non_positive_roots_1)
    
    return (yan_val_storer_1,empty_HiC_list,r_1_list_1,non_positive_roots_1)

if __name__ == '__main__':
    list_trial,list_trials = Region_partition(0,0.6,200)
    spaced = np.linspace(0,0.6,200)
    print(list_trial, spaced)#list_trials)
'''