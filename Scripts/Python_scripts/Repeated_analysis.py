#First let us load in our HiC matrix:
import numpy as np
from iced import normalization
import csv

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

def submatrix_reader(matrix,submatrix_shape,stride = 1,offset = (0,0),pattern = "diagonal"):
    #Logic should be working only for symmetric submatrices, not properly generalized to non square matrices
    submatrix_no = min(matrix.shape) - max(submatrix_shape) + 1 - max(offset)
    submatrices=[]
    for i in range(0,submatrix_no,stride):
        submatrix = matrix[i+offset[0]:i+submatrix_shape[0]+offset[0],i+offset[1]:i+submatrix_shape[1] + offset[1]]
        submatrices.append(submatrix)
    return submatrices

#https://stackoverflow.com/questions/38708621/how-to-calculate-percentage-of-sparsity-for-a-numpy-array-matrix
def sparsity_calc(matrix):
    matrix = np.nan_to_num(matrix,0)
    sparsity = 1.0 - ( np.count_nonzero(matrix) / float(matrix.size) )
    return sparsity

HiC_data = matrix_reader("C:/Users/HP/Downloads/Siddharth/Data/hg38/GM12878/100kb/raw/chr1.mat.txt")
np.fill_diagonal(HiC_data,0) #Making diagonals of our matrix zero
Regions_to_check = (143000001,158000001)
Resolution = 100000
bins_to_check = (Regions_to_check[0]//Resolution,Regions_to_check[1]//Resolution)
HiC_data_to_work = HiC_data[bins_to_check[0]:bins_to_check[1],bins_to_check[0]:bins_to_check[1]]
ICEd_mat = normalization.ICE_normalization(HiC_data_to_work)
logged_ICE_mat = np.log(ICEd_mat + 1)
logged_ICEd_submats = submatrix_reader(logged_ICE_mat,(10,10))
submat_sum_list = list(map(np.sum,logged_ICEd_submats))
submat_sparsity = list(map(sparsity_calc,logged_ICEd_submats))
print(submat_sum_list,'\n',submat_sparsity)

K_50_ICE_path = "C:/Users/HP/Downloads/Siddharth/Scripts/Python_scripts/GM12878_ICE_K_50_trial_list.pickle"