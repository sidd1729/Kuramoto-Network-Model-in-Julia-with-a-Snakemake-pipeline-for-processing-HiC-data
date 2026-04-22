import numpy as np
import scipy
import csv
import subprocess
from multiprocessing import Pool
import os
from iced import normalization
import matplotlib.pyplot as plt
from tqdm import tqdm

def submatrix_reader(matrix,submatrix_shape,stride = 1,offset = (0,0),pattern = "diagonal"):
    #Logic should be working only for symmetric submatrices, not properly generalized to non square matrices
    submatrix_no = min(matrix.shape) - max(submatrix_shape) + 1 - max(offset)
    submatrices=[]
    for i in range(0,submatrix_no,stride):
        submatrix = matrix[i+offset[0]:i+submatrix_shape[0]+offset[0],i+offset[1]:i+submatrix_shape[1] + offset[1]]
        submatrices.append(submatrix)
    return submatrices

def matrix_reader(filename):
    file = open('{}'.format(filename), 'r')
    content = file.readlines()
    #print(len(content[1]))
    file.close()

    tsv_reader = csv.reader(content, delimiter='\t')
    counter = 0
    row_list = []
    for row in tsv_reader :
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

def list_writer(filename,variable_name):
    file = open('{}'.format(filename), 'w')
    for element in variable_name:
            file.write(str(element) +'\t')
    file.close()

def tuple_list_writer(filename,variable_name):
    file = open('{}'.format(filename), 'w')
    for t in variable_name:
      file.write(' '.join(str(s) for s in t) + '\n')
    file.close()

def K_50_val_calc(HiC_data_file,Regions_to_check,Resolution,node_no,data_root,Matrix_Form,Order_parameter_form,to_write = False):
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

def simple_plotter(yan_val_storer,image_root,Matrix_Form,Order_parameter_form,Chromosome_no,Resolution,to_write = False):
    temp_1 = []
    temp_2 = []
    for i,j in yan_val_storer:
        temp_1.append(j)
        temp_2.append(i)
    plt.plot(temp_1,temp_2,'o',markersize = 0.9,color = "blue",label = "{} Kuramoto Order Parameter".format(Order_parameter_form))
    plt.legend()
    plt.title("K_0.5 value ({}) vs chromosome {} position for {} size bins".format(Matrix_Form,Chromosome_no,Resolution))
    plt.xlabel("Chromosome {} Position for {} size bin".format(Chromosome_no,Resolution))
    plt.ylabel("K value for r = 0.5")
    Img_Dir_exist = os.path.exists("{}/{}/{}".format(image_root,Matrix_Form,Order_parameter_form))
    if not Img_Dir_exist:
         os.makedirs("{}/{}/{}".format(image_root,Matrix_Form,Order_parameter_form))
    if to_write == True:
        plt.savefig("{}/{}/{}/{}_K_0_5".format(image_root,Matrix_Form,Order_parameter_form,Order_parameter_form))
    plt.close()

def K_50_calc_plotter(info):
    HiC_data_file = info[0]
    #np.fill_diagonal(HiC_data,0) #Making diagonals of our matrix zero
    Regions_to_check = info[1]
    Resolution = info[2]
    node_no = info[3]
    data_root = info[4]
    image_root = info[5]
    Chromosome_no = info[6]
    Matrix_Form = info[7]
    Order_parameter_form = info[8]
   #yan_val_storer_1,empty_HiC_list,r_1_list_1,non_positive_roots_1  =  K_50_val_calc(HiC_data_file,Regions_to_check,Resolution,node_no,data_root,Matrix_Form,Order_parameter_form)
    yan_val_storer_1,empty_HiC_list,r_1_list_1,non_positive_roots_1  =  K_50_val_calc_R(HiC_data_file,Regions_to_check,Resolution,node_no,data_root,Matrix_Form,Order_parameter_form)
    simple_plotter(yan_val_storer_1,image_root,Matrix_Form,Order_parameter_form,Chromosome_no,Resolution)
    
def K_50_val_calc_R(HiC_data_file,Regions_to_check,Resolution,node_no,data_root,Matrix_Form,Order_parameter_form,to_write = False):
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
    if Order_parameter_form == "Kuramoto":
        k_super_storer_path =   "{}/{}/{}_k_values_kura_order".format(data_root,Matrix_Form,Matrix_Form)
        r_super_storer_path  =  "{}/{}/{}_r_values_kura_order".format(data_root,Matrix_Form,Matrix_Form)
    elif Order_parameter_form == "Universal":
        k_super_storer_path =  "{}/{}/{}_k_values_uni_order".format(data_root,Matrix_Form,Matrix_Form)
        r_super_storer_path =  "{}/{}/{}_r_values_uni_order".format(data_root,Matrix_Form,Matrix_Form)
    yan_val_storer_1 = []
    empty_HiC_list =[]
    r_1_list_1 = []
    non_positive_roots_1 = []
    results= subprocess.check_output(['RScript',"C:/Users/HP/Downloads/Siddharth/Scripts/R_scripts/Curve_pipeline.R","{}".format(k_super_storer_path ),"{}".format(r_super_storer_path)],text=True)
    results_list = results.split("\n")
    del(results_list[len(results_list)-1])
    final_output = tuple(map(lambda x : float(x[4:]), results_list)) #Output in the form of K_50 value,R_value
    filtered_root = final_output[0]
    i = final_output[3]
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
    info_list = []
    Matrix_Form_List = ["Raw","ICE","Z"]
    Order_Parameter_Form_List = ["Kuramoto","Universal"]
    info_list_root = [["C:/Users/HP/Downloads/Siddharth/Data/hg38/GM12878/100kb/raw/chr1.mat.txt",(143000001,158000001),100000,10,"C:/Users/HP/Downloads/Siddharth/Data/hg38/GM12878/100kb/Chromosome_1","C:/Users/HP/Downloads/Siddharth/Thesis_Figures/Weighted_Network/GM12878/100kb/Chromosome_1",1],["C:/Users/HP/Downloads/Siddharth/Data/hg38/GM12878/10kb/raw/chr1.mat.txt",(143000001,158000001),10000,100,"C:/Users/HP/Downloads/Siddharth/Data/hg38/GM12878/10kb/Chromosome_1","C:/Users/HP/Downloads/Siddharth/Thesis_Figures/Weighted_Network/GM12878/10kb/Chromosome_1",1],
    ["C:/Users/HP/Downloads/Siddharth/Data/hg38/stem_cell_H1_HESC/100kb/raw/chr1.mat.txt",(143000001,158000001),100000,10,"C:/Users/HP/Downloads/Siddharth/Data/hg38/stem_cell_H1_HESC/100kb/Chromosome_1","C:/Users/HP/Downloads/Siddharth/Thesis_Figures/Weighted_Network/stem_cell_H1_HESC/100kb/Chromosome_1",1]]
    for Matrix_Form in Matrix_Form_List:
        for Order_parameter_form in Order_Parameter_Form_List:
            for info_list_root_segment in info_list_root:
                #print(Matrix_Form,Order_parameter_form)
                info_list_root_segment_copy = info_list_root_segment.copy()
                #print(info_list_root_segment_copy)
                info_list_root_segment_copy.extend([Matrix_Form,Order_parameter_form])
                #print(info_list_root_segment)
                info_list.append(tuple(info_list_root_segment_copy))
                #print(info_list)

    for i in info_list:
        print(i[7:9])
    #info_list = [("C:/Users/HP/Downloads/Siddharth/Data/hg38/GM12878/100kb/raw/chr1.mat.txt",(143000001,158000001),100000,10,"C:/Users/HP/Downloads/Siddharth/Data/hg38/GM12878/100kb/Chromosome_1","C:/Users/HP/Downloads/Siddharth/Thesis_Figures/Weighted_Network/GM12878/100kb/Chromosome_1",1,"Raw","Kuramoto"),("C:/Users/HP/Downloads/Siddharth/Data/hg38/GM12878/100kb/raw/chr1.mat.txt",(143000001,158000001),100000,10,"C:/Users/HP/Downloads/Siddharth/Data/hg38/GM12878/100kb/Chromosome_1","C:/Users/HP/Downloads/Siddharth/Thesis_Figures/Weighted_Network/GM12878/100kb/Chromosome_1",1,"Raw","Universal")]
    processes_to_use = 1
    max_ = len(info_list)
    print(max_)
    #print(info_list)
    
    with Pool(processes=processes_to_use) as p, tqdm(total=max_) as pbar:
            for result in p.imap(K_50_calc_plotter, info_list):
                pbar.update()
                pbar.refresh()

            
        