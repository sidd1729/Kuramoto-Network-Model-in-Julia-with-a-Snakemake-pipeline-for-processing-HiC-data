import sys #Reading in command line arguments
import csv #Readomg in tab separated files
import numpy as np #Matrix/Array operations
import pyarrow.feather as feather #Reading and Writing data to format which can be used in both python and R
import matplotlib.pyplot as plt #Plotting graphs
from iced import normalization #ICED normalization of HiC data
import cmasher #Better perceptually uniform colormaps to visualize HiC data
import re #Regular expressions for parsing text files
import pandas as pd #Easier processing of parsed text files using pandas dataframes
from sklearn.decomposition import PCA #Performing principal component analysis of HiC matrices
from tqdm import tqdm #Adding a timer to for loops/other operations on iterables

#Function to read in the HiC matrices (from text/npy formats)
def matrix_reader(filename,delete = False):
    if filename[-3:] == "npy":
        return np.load(filename)
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

#Function to read submatrices of a given shape from a matrix:
def submatrix_reader(matrix,submatrix_shape,stride = 1,offset = (0,0),pattern = "diagonal"):
    #Logic should be working only for symmetric submatrices, not properly generalized to non square matricess
    submatrix_no = min(matrix.shape) - max(submatrix_shape) + 1 - max(offset)
    submatrices=[]
    for i in range(0,submatrix_no,stride):
        submatrix = matrix[i+offset[0]:i+submatrix_shape[0]+offset[0],i+offset[1]:i+submatrix_shape[1] + offset[1]]
        submatrices.append(submatrix)
    return submatrices

#Function to apply offset to a list and converts it into a numpy array for ease of use:
def offset_applier(original_list,offset):
     return np.array(list(map(lambda x: x+offset, original_list)))

#Function to subtract start bin number for plotting purposes:
def start_bin_sub(original_list, start_bins_to_check):
     return np.array(list(map(lambda x: x-start_bins_to_check,original_list)))

#Function to count the number of characters in a fasta file(Relevant for gene density calculations):
def file_char_counter(filename):
    file = open(filename,'r')
    content = file.read()
    file.close()
    len_N = len(re.findall(r'[N]',content))
    len_char = len(re.findall(r'[a-zA-Z]',content))
    return (len_N,len_char)

#Function to calculate insulation scores from a HiC matrix
def insulation_score_calc(HiC_matrix,node_no=5,stride = 1,offset = (0,5)):
        insulated_sub =  submatrix_reader(HiC_matrix,(node_no,node_no),stride,offset)
        insulation_score = np.array(list(map(np.sum,insulated_sub)))
        return insulation_score

#Function to calculate Log 2 local - distal interaction ratio:
#Obtaining local distal interaction ratio:
def local_distal_interaction_ratio(HiC_matrix,Resolution):
        Distal_Interaction_List = [0]*HiC_matrix.shape[0]
        Local_Interaction_List = [0]*HiC_matrix.shape[0]

        for i in range(0,HiC_matrix.shape[0]):
            for j in range(0,HiC_matrix.shape[1]):
                if abs(j-i) < 3000000/Resolution :
                    Local_Interaction_List[i]+=HiC_matrix[i,j]
                else:
                    Distal_Interaction_List[i]+=HiC_matrix[i,j]
        log_2_LDR_list = [0]*HiC_matrix.shape[0]
        for i in range(0,HiC_matrix.shape[0]):
            if Distal_Interaction_List[i] == 0 :
                Distal_Interaction_List[i] += 0.0000001
            log_2_LDR_list[i] = np.log2(Local_Interaction_List[i]/(Distal_Interaction_List[i]))
        return log_2_LDR_list

#Function for calculating Gene Density
def RefSeq_Gene_Density_Calc(filename,bin_size,number_characters_in_fasta_file):
    file = open(filename,'r')
    content = file.readlines()
    file.close()
    tsv_reader = csv.reader(content, delimiter='\t')
    counter = 0
    row_list = []
    for row in tsv_reader :
            #del row[len(row)-1]
            for i in range(0,len(row)):
                 row[i] = str(row[i])
            counter += 1
            row_list.append(row)
    j = pd.DataFrame(row_list)
    j.columns = j.iloc[0,:] 
    j.drop(0)
    j['txStart'] = pd.to_numeric(j['txStart'], downcast='integer', errors='coerce')
    j['txEnd'] = pd.to_numeric(j['txEnd'], downcast='integer', errors='coerce')
    bin_list = []
    gene_density_list = []
    for i in range(0,number_characters_in_fasta_file,bin_size):
        bin_list.append((i,i+bin_size))
        #filtered = j.loc[((j['txStart'] > i) & (j['txEnd'] < i+bin_size)) | ((j['txStart'] > i) & (j['txStart'] < i+bin_size) & (j['txEnd'] < i+2*bin_size))]# Select entries lying in a bin, current ignoring those spanning between bins
        filtered = j.loc[(j['txStart'] > i) & (j['txStart'] < i+bin_size)]
        unique_name_in_filtered =  set(filtered.loc[:,'name2'])
        #print(filtered,i,i+bin_size)
        gene_density_list.append(len(unique_name_in_filtered))

    #plot_list = [i for i in range(0,len(gene_density_list))]
    return (gene_density_list)

#Function to read wig files of ChiP seq marks, obtained after downloading bam files from ENCODE, generating index/bai files using Rsamtools and using bamCoverage from deeptools
def wig_reader(filename):
    chr_no_list = []
    start_no_list = []
    stop_no_list = []
    count_list = []
    with open(filename, 'r') as file:
        # Reading and printing the entire file line by line
        while True:
            line = file.readline()
            if re.match(r'#',line):
                continue
            #else:
                #print(line,line.split('\t'))
                #print(re.match(r'#',line))
            elif re.match(r'([a-z]{3}[0-9{1-2}|X|Y])\t([0-9]+)\t([0-9]+)\t([0-9]+)',line):
                result = re.match(r'([a-z]{3}[0-9{1-2}|X|Y])\t([0-9]+)\t([0-9]+)\t([0-9]+)',line)
                print(result.groups())
                result_groups = result.groups()
                #if (result_groups[0] == "{}".format(desired_chr)) and (int(result_groups[1]) >= desired_region[0]) and (int(result_groups[2]) <= desired_region[1]):
                chr_no_list.append(result_groups[0])
                start_no_list.append(result_groups[1])
                stop_no_list.append(result_groups[2])
                count_list.append(result_groups[3])
            if not line:
                break
    start_no_list = list(map(int,start_no_list))
    stop_no_list = list(map(int,stop_no_list))
    count_list = list(map(int,count_list))
    return (start_no_list,stop_no_list,count_list)

#Function pad 0 in ChiP sequencing result lists:
def ChiP_seq_padder(start_no_list,count_list,Resolution = 100000,chromosome_limit = 248956425):
    #value_list = []
    for value in range(0,248956425,100000):
             if value not in start_no_list:
                    count_list.insert((value//Resolution),0)
                    #value_list.append(value)
    #If you wish to verify that the program is working correctly then (uncomment all the lines involving value list and you should see that the below print only gives 0 values)
    #print(np.array(count_list)[np.array(list(map(lambda x : x//100000,values_list)))])
    return count_list
#Function to read in mappability scores/bedgraphs downloaded from https://bismap.hoffmanlab.org/ after they have been unzipped:
def bedgraph_reader(filename,explicit_chrom = None,binned = False,bin_size = None,chr_size = None):
    #filename = r"\\wsl.localhost\Ubuntu\home\default\deeptools\HESC_H3K27.wig"
    chr_no_list = []
    start_no_list = []
    stop_no_list = []
    count_list = []
    #desired_chr = "chr1"
    #desired_region = (158000001,173000001)
    with open(filename, 'r') as file:
        # Reading and printing the entire file line by line
        while True:
            line = file.readline()
            if re.match(r'#',line):
                continue
            #else:
                #print(line,line.split('\t'))
                #print(re.match(r'#',line))
            elif re.match(r'([a-z]{3}[0-9{1-2}|X|Y])\t([0-9]+)\t([0-9]+)\t([0|1]\.(?:[0-9]+))',line):
                if explicit_chrom == None:
                    result = re.match(r'([a-z]{3}[0-9{1-2}|X|Y])\t([0-9]+)\t([0-9]+)\t([0|1]\.(?:[0-9]+))',line)
                else:
                    result = re.match(f'({explicit_chrom})\t([0-9]+)\t([0-9]+)\t([0|1]\\.(?:[0-9]+))',line)
                if result == None:
                    break
                #print(result.groups())
                result_groups = result.groups()
                #if (result_groups[0] == "{}".format(desired_chr)) and (int(result_groups[1]) >= desired_region[0]) and (int(result_groups[2]) <= desired_region[1]):
                chr_no_list.append(result_groups[0])
                start_no_list.append(result_groups[1])
                stop_no_list.append(result_groups[2])
                count_list.append(result_groups[3])
                
            if not line:
                break
    start_no_list = list(map(int,start_no_list))
    stop_no_list = list(map(int,stop_no_list))
    count_list = list(map(float,count_list))
    if binned == False:
        processing_df = pd.DataFrame({'Start':start_no_list,'Stop':stop_no_list,'Score':count_list})
        print(processing_df)
        bin_list=[]
        count_list = []
        for i in tqdm(range(0,chr_size,bin_size)):
            bin_list.append((i,i+bin_size))
            filtered = processing_df.loc[(processing_df['Start'] >= i) & (processing_df['Start'] <= i+bin_size)]
            count_list.append(filtered['Score'].sum())
        #print(filtered,i,i+bin_size)
        #count_list.append(len(unique_name_in_filtered))
        start_no_list = list(map(lambda x: x[0],bin_list))
        stop_no_list = list(map(lambda x: x [1],bin_list))
    return (start_no_list,stop_no_list,count_list)

#Relevant parameters:
window_size = 5
offset = np.floor((window_size+1)/2)
try:
    Regions_to_check = (int(sys.argv[1]),int(sys.argv[2]))
except:
    Regions_to_check = (143000001,158000001)
Resolution = 100000
bins_to_check = (Regions_to_check[0]//Resolution,Regions_to_check[1]//Resolution)
Entire_range = np.arange(0,int(2490),1)

#Reading in the relevant matrix and applying (Log (x+1) transformation(ICE_ normalization(matrix))
HiC_data = matrix_reader("/home/ksslab/Siddharth/Program_Directory/Data/hg38/GM12878/100kb/raw/chr1.mat.txt")
np.fill_diagonal(HiC_data,0) #Making diagonals of our matrix zero
HiC_data_to_work = HiC_data[bins_to_check[0]:bins_to_check[1]+1,bins_to_check[0]:bins_to_check[1]+1]
ICEd_data = normalization.ICE_normalization(HiC_data)
ICEd_data_mat = ICEd_data[bins_to_check[0]:bins_to_check[1]+1,bins_to_check[0]:bins_to_check[1]+1] #For LDR calculation
logged_overall_ICE_mat = np.log(ICEd_data+1)
logged_ICE_mat = logged_overall_ICE_mat[bins_to_check[0]:bins_to_check[1]+1,bins_to_check[0]:bins_to_check[1]+1]


#Repeating the same for HESC data
HiC_data_HESC = matrix_reader("/home/ksslab/Siddharth/Program_Directory/Data/hg38/stem_cell_H1_HESC/100kb/raw/chr1.mat.txt")
np.fill_diagonal(HiC_data_HESC,0) #Making diagonals of our matrix zero
HiC_data_to_work_HESC = HiC_data_HESC[bins_to_check[0]:bins_to_check[1]+1,bins_to_check[0]:bins_to_check[1]+1]
ICEd_data_HESC = normalization.ICE_normalization(HiC_data_HESC)
ICEd_data_mat_HESC = ICEd_data_HESC[bins_to_check[0]:bins_to_check[1]+1,bins_to_check[0]:bins_to_check[1]+1]
logged_overall_ICE_mat_HESC = np.log(ICEd_data_HESC+1)
logged_ICE_mat_HESC = logged_overall_ICE_mat_HESC[bins_to_check[0]:bins_to_check[1]+1,bins_to_check[0]:bins_to_check[1]+1]



#Perpetually 0 indices to be ignored,read in and corrected:
perpet_1_indice_arr = (feather.read_feather("/home/ksslab/Siddharth/Program_Directory/Data/Simulation_Results/Linux/GM12878/Chr1_perpet_1_Kura.feather")).to_numpy()
perpet_1_vals = np.array(list(map(lambda x : x[0][0],perpet_1_indice_arr)))
problem_indices_cell_type_1 = perpet_1_vals[perpet_1_vals != -1]
problem_indices_cell_type_1 = np.array(list(map(lambda x: x-1,problem_indices_cell_type_1)))
Corrected_range = np.setdiff1d(Entire_range, problem_indices_cell_type_1)
Truly_corrected_range = np.intersect1d(Corrected_range[bins_to_check[0]-1<Corrected_range],Corrected_range[Corrected_range < bins_to_check[1]-window_size+1])
Truly_corrected_range = start_bin_sub(Truly_corrected_range,bins_to_check[0])
#print(Truly_corrected_range)

#Perpetually 0 indices to be ignored,read in and corrected(done for cell_type 2)
perpet_1_indice_arr_2 = (feather.read_feather("/home/ksslab/Siddharth/Program_Directory/Data/Simulation_Results/Linux/stem_cell_H1_HESC/Chr1_perpet_1_Kura.feather")).to_numpy()
perpet_1_vals_2 = np.array(list(map(lambda x : x[0][0],perpet_1_indice_arr_2)))
problem_indices_cell_type_2 = perpet_1_vals_2[perpet_1_vals_2 != -1]
problem_indices_cell_type_2 = np.array(list(map(lambda x: x-1,problem_indices_cell_type_2)))
Corrected_range_2 = np.setdiff1d(Entire_range, problem_indices_cell_type_2)
Truly_corrected_range_2 = np.intersect1d(Corrected_range_2[bins_to_check[0]-1<Corrected_range_2],Corrected_range_2[Corrected_range_2 < bins_to_check[1]-window_size+1])
Truly_corrected_range_2 = start_bin_sub(Truly_corrected_range_2,bins_to_check[0])

#List with indices to do plotting:
appropriate_range_plotting_list = list(range(bins_to_check[0],bins_to_check[1]+1))
appropriate_range_plotting_list = start_bin_sub(appropriate_range_plotting_list,bins_to_check[0])
#print(len(appropriate_range_plotting_list),'\n',appropriate_range_plotting_list)

#Making offset versions of previously generated lists:
appropriate_range_plotting_list_offset  =offset_applier(appropriate_range_plotting_list,offset)
Truly_corrected_range_offset = offset_applier(Truly_corrected_range,offset)
Truly_corrected_range_offset_2 = offset_applier(Truly_corrected_range_2,offset)

#print(appropriate_range_plotting_list,'\n',appropriate_range_plotting_list_offset,'\n',Truly_corrected_range,'\n',Truly_corrected_range_offset)

#Reading in slice of plucked array:
plucked = feather.read_feather("/home/ksslab/Siddharth/Program_Directory/Data/Simulation_Results/Linux/GM12878/Chr1_K_50_Kura.feather")
plucked_arr = plucked.to_numpy()
Whole_chr_1_K_50 = np.array(list(map(lambda x : x[0][0][0],plucked_arr)))

 #Performing the reading of another cell type:
plucked_HESC = feather.read_feather("/home/ksslab/Siddharth/Program_Directory/Data/Simulation_Results/Linux/stem_cell_H1_HESC/Chr1_K_50_Kura.feather")
plucked_arr_HESC = plucked_HESC.to_numpy()
Whole_chr_1_K_50_HESC = np.array(list(map(lambda x : x[0][0][0],plucked_arr_HESC)))


#Other parameters/Calculations to be read in :

#Insulation score:
insulated_ICEd_score = insulation_score_calc(logged_ICE_mat)
insulated_ICEd_score_2 = insulation_score_calc(logged_ICE_mat_HESC)
#GM12878 Mat:

#LDR:
LDR_score = local_distal_interaction_ratio(ICEd_data_mat,Resolution)
LDR_score_HESC = local_distal_interaction_ratio(ICEd_data_mat_HESC,Resolution)

#Gene Density (Cell type independent)
    #First counting number of characters in corresponding fasta file:
len_N,len_char = file_char_counter("/home/ksslab/Siddharth/Program_Directory/Data/hg38/chr1.fa")
Gene_density_calc = RefSeq_Gene_Density_Calc("/home/ksslab/Siddharth/Program_Directory/Data/hg38/Curated_ref_seq_chr1.txt",100000,len_char)

#H3K27Me3
GM12878_chr1_K27Me_count_list = wig_reader("/home/ksslab/Siddharth/Program_Directory/Data/ChiP_Seq_Data/GM12878/H3K27Me3/GM12878_Chr1_H3K27Me3.wig")
GM12878_chr1_K27Me_count_list = ChiP_seq_padder(GM12878_chr1_K27Me_count_list[0],GM12878_chr1_K27Me_count_list[2])
HESC_chr1_K27Me_count_list = wig_reader("/home/ksslab/Siddharth/Program_Directory/Data/ChiP_Seq_Data/H1_HESC/H3K27Me3/HESC_Chr1_H3K27Me3.wig")
HESC_chr1_K27Me_count_list = ChiP_seq_padder(HESC_chr1_K27Me_count_list[0],HESC_chr1_K27Me_count_list[2])


#H3K4Me3
GM12878_chr1_K4Me_count_list = wig_reader("/home/ksslab/Siddharth/Program_Directory/Data/ChiP_Seq_Data/GM12878/H3K4Me3/GM12878_Chr1_H3K4Me3.wig")
GM12878_chr1_K4Me_count_list = ChiP_seq_padder(GM12878_chr1_K4Me_count_list[0],GM12878_chr1_K4Me_count_list[2])
HESC_chr1_K4Me_count_list = wig_reader("/home/ksslab/Siddharth/Program_Directory/Data/ChiP_Seq_Data/H1_HESC/H3K4Me3/HESC_Chr1_H3K4Me3.wig")
HESC_chr1_K4Me_count_list = ChiP_seq_padder(HESC_chr1_K4Me_count_list[0],HESC_chr1_K4Me_count_list[2])

#H3K27Ac
GM12878_chr1_K27Ac_count_list = wig_reader("/home/ksslab/Siddharth/Program_Directory/Data/ChiP_Seq_Data/GM12878/H3K27Ac/GM12878_Chr1_H3K27Ac.wig")
GM12878_chr1_K27Ac_count_list = ChiP_seq_padder(GM12878_chr1_K27Ac_count_list[0],GM12878_chr1_K27Ac_count_list[2])
HESC_chr1_K27Ac_count_list = wig_reader("/home/ksslab/Siddharth/Program_Directory/Data/ChiP_Seq_Data/H1_HESC/H3K27Ac/HESC_Chr_1_H3K27Ac.wig")
HESC_chr1_K27Ac_count_list = ChiP_seq_padder(HESC_chr1_K27Ac_count_list[0],HESC_chr1_K27Ac_count_list[2])

#CTCF
GM12878_chr1_CTCF_count_list = wig_reader("/home/ksslab/Siddharth/Program_Directory/Data/ChiP_Seq_Data/GM12878/CTCF/GM12878_Chr1_CTCF.wig")
GM12878_chr1_CTCF_count_list = ChiP_seq_padder(GM12878_chr1_CTCF_count_list[0],GM12878_chr1_CTCF_count_list[2])
HESC_chr1_CTCF_count_list = wig_reader("/home/ksslab/Siddharth/Program_Directory/Data/ChiP_Seq_Data/H1_HESC/CTCF/HESC_Chr1_CTCF.wig")
HESC_chr1_CTCF_count_list = ChiP_seq_padder(HESC_chr1_CTCF_count_list[0],HESC_chr1_CTCF_count_list[2])


#Chromosome 1, all the ChiP seq readouts

#PCs
no_components = 3
pca = PCA(n_components = no_components)
PCA_output = pca.fit(logged_overall_ICE_mat)
ICEd_PCs = PCA_output.components_

PCA_output_HESC = pca.fit(logged_overall_ICE_mat_HESC)
ICEd_PCs_HESC = PCA_output_HESC.components_

#PC1
GM_12878_PC1 = ICEd_PCs[0,:][bins_to_check[0]:bins_to_check[1]+1]
HESC_PC1 = ICEd_PCs_HESC[0,:][bins_to_check[0]:bins_to_check[1]+1]

#PC2
GM_12878_PC2 = ICEd_PCs[1,:][bins_to_check[0]:bins_to_check[1]+1]
HESC_PC2 = ICEd_PCs_HESC[1,:][bins_to_check[0]:bins_to_check[1]+1]

#PC3
GM_12878_PC3 = ICEd_PCs[2,:][bins_to_check[0]:bins_to_check[1]+1]
HESC_PC3 = ICEd_PCs_HESC[2,:][bins_to_check[0]:bins_to_check[1]+1]

#Mappability (Genome property)
Mappability_count_list = bedgraph_reader("/home/ksslab/Siddharth/Program_Directory/Data/hg38/mmap_k100.umap.bedgraph","chr1",False,100000,248956425)[2]

fig = plt.figure(figsize=(20,20))
gs = fig.add_gridspec(13,2, hspace=0,height_ratios = [2,10,1,1,1,1,1,1,1,1,1,1,1],width_ratios =[1,1])
axs = gs.subplots(sharex = True)

#K_50 values without r = 1
axs[0,0].plot(Truly_corrected_range_offset,Whole_chr_1_K_50[Truly_corrected_range],'o-',markersize = 4,color = 'orange',label = "K_50 without r = 1")
axs[0,0].legend(fontsize = 6)

#Log transformed ICE normalized HiC matrixiC_data_to_work_HESC
axs[1,0].imshow(logged_ICE_mat,cmap = 'cmr.sunburst_r',aspect = "auto")
#axs[1,]

#Insulation score
axs[2,0].plot(appropriate_range_plotting_list_offset[0:len(insulated_ICEd_score)],insulated_ICEd_score,color = 'orange',label = "Insulation score")
axs[2,0].legend(fontsize = 6)

#Log 2 Local - Distal Interaction Ratio
axs[3,0].plot(appropriate_range_plotting_list,LDR_score,color = 'orange',label = "Log 2 Local Distal Interaction Ratio")
axs[3,0].legend(fontsize = 6)
#plt.plot()

#Gene Density
axs[4,0].plot(appropriate_range_plotting_list,Gene_density_calc[bins_to_check[0]:bins_to_check[1]+1],color = 'orange',label = "Gene Density")
axs[4,0].legend(fontsize = 6)

#H3K27Me3
axs[5,0].plot(appropriate_range_plotting_list,GM12878_chr1_K27Me_count_list[bins_to_check[0]:bins_to_check[1]+1],color = 'orange', label = "GM12878 H3K27Me3")
axs[5,0].legend(fontsize = 6)

#H3K27Ac
axs[6,0].plot(appropriate_range_plotting_list,GM12878_chr1_K27Ac_count_list[bins_to_check[0]:bins_to_check[1]+1],color = 'orange', label = "GM12878 H3K27Ac")
axs[6,0].legend(fontsize = 6)

#H3K4Me3
axs[7,0].plot(appropriate_range_plotting_list,GM12878_chr1_K4Me_count_list[bins_to_check[0]:bins_to_check[1]+1],color = 'orange', label = "GM12878 H3K4Me3")
axs[7,0].legend(fontsize = 6)

#CTCF
axs[8,0].plot(appropriate_range_plotting_list,GM12878_chr1_CTCF_count_list[bins_to_check[0]:bins_to_check[1]+1],color = 'orange',label = "GM12878 CTCF")
axs[8,0].legend(fontsize = 6)

#PC1
axs[9,0].plot(appropriate_range_plotting_list,GM_12878_PC1,color = 'orange',label = "GM12878 PC1")
axs[9,0].legend(fontsize = 6)

#PC2
axs[10,0].plot(appropriate_range_plotting_list,GM_12878_PC2,color = 'orange',label = "GM12878 PC2")
axs[10,0].legend(fontsize = 6)

#PC3
axs[11,0].plot(appropriate_range_plotting_list,GM_12878_PC3,color = 'orange',label = "GM12878 PC3")
axs[11,0].legend(fontsize = 6)

#Mappability
axs[12,0].plot(appropriate_range_plotting_list,Mappability_count_list[bins_to_check[0]:bins_to_check[1]+1],color = 'orange', label = "Mappability scores")
axs[12,0].legend(fontsize = 6)

#K_50 values without r = 1
axs[0,1].plot(Truly_corrected_range_offset_2,Whole_chr_1_K_50_HESC[Truly_corrected_range_2],'o-',markersize = 4,color = 'orange',label = "K_50 HESC without r = 1")
axs[0,1].legend(fontsize = 6)

#Log transformed ICE normalized HiC matrix
axs[1,1].imshow(logged_ICE_mat_HESC,cmap = 'cmr.sunburst_r',aspect = "auto")
#axs[1,]

#Insulation score
axs[2,1].plot(appropriate_range_plotting_list_offset[0:len(insulated_ICEd_score_2)],insulated_ICEd_score_2,color = 'orange',label = "Insulation score HESC")
axs[2,1].legend(fontsize = 6)

#Log 2 Local - Distal Interaction Ratio
axs[3,1].plot(appropriate_range_plotting_list,LDR_score_HESC,color = 'orange',label = "Log 2 Local Distal Interaction Ratio HESC")
axs[3,1].legend(fontsize = 6)
#plt.plot()

#Gene Density
axs[4,1].plot(appropriate_range_plotting_list,Gene_density_calc[bins_to_check[0]:bins_to_check[1]+1],color = 'orange',label = "Gene Density")
axs[4,1].legend(fontsize = 6)

#H3K27Me3
axs[5,1].plot(appropriate_range_plotting_list,HESC_chr1_K27Ac_count_list[bins_to_check[0]:bins_to_check[1]+1],color = 'orange', label = "HESC H3K27Me3")
axs[5,1].legend(fontsize = 6)

#H3K27Ac
axs[6,1].plot(appropriate_range_plotting_list,HESC_chr1_K27Ac_count_list[bins_to_check[0]:bins_to_check[1]+1],color = 'orange', label = "HESC H3K27Ac")
axs[6,1].legend(fontsize = 6)

#H3K4Me3
axs[7,1].plot(appropriate_range_plotting_list,HESC_chr1_K4Me_count_list[bins_to_check[0]:bins_to_check[1]+1],color = 'orange', label = "HESC H3K4Me3")
axs[7,1].legend(fontsize = 6)

#CTCF
axs[8,1].plot(appropriate_range_plotting_list,HESC_chr1_CTCF_count_list[bins_to_check[0]:bins_to_check[1]+1],color = 'orange',label = "HESC CTCF")
axs[8,1].legend(fontsize = 6)

#PC1
axs[9,1].plot(appropriate_range_plotting_list,HESC_PC1,color = 'orange',label = "HESC PC1")
axs[9,1].legend(fontsize = 6)

#PC2
axs[10,1].plot(appropriate_range_plotting_list,HESC_PC2,color = 'orange',label = "HESC PC2")
axs[10,1].legend(fontsize = 6)

#PC3
axs[11,1].plot(appropriate_range_plotting_list,HESC_PC3,color = 'orange',label = "HESC PC3")
axs[11,1].legend(fontsize = 6)


#Mappability
axs[12,1].plot(appropriate_range_plotting_list,Mappability_count_list[bins_to_check[0]:bins_to_check[1]+1],color = 'orange', label = "Mappability scores")
axs[12,1].legend(fontsize = 6)

#plt.savefig(f"/home/ksslab/Siddharth/Temp_plot_dump/K_50_params_plot_{int(sys.argv[1])}_{int(sys.argv[2])}.png")
plt.show()
plt.close()