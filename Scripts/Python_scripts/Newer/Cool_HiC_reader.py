import cooler
import cooltools
import matplotlib.pyplot as plt
import hicstraw
import numpy as np
import pandas as pd
import csv
import pyarrow.feather as feather
import itertools 
import os
import h5py
import bioframe
#help(cooler.Cooler)

clr = cooler.Cooler("/home/ksslab/Siddharth/Program_Directory/Data/Actual_Data/Raw_cool/40kb/12_hpi_S165.cool")
print(f'chromosomes: {clr.chromnames}, binsize: {clr.binsize}')

chromstarts = []
for i in clr.chromnames:
    print(f'{i} : {clr.extent(i)}')
    chromstarts.append(clr.extent(i)[0])

np_mat = clr.matrix(balance=False).fetch(('1',0,59_578_282)) #15578282))#
print(type(np_mat))
plt.imshow(np_mat)
#plt.savefig("/home/ksslab/Siddharth/Temp_plot_dump/Cool_tools_chr_1_0_15178282.png")
plt.show()
plt.close()

complete_cooler_mat = clr.matrix(balance=False)

hic = hicstraw.HiCFile("/home/ksslab/Siddharth/Program_Directory/Data/Actual_Data/Raw_hic/40kb/12_hpi_S165.hic")
print(hic.getChromosomes())
print(hic.getGenomeID())
print(hic.getResolutions())

#help(hic.getMatrixZoomData)
mzd = hic.getMatrixZoomData('1', '1', "observed", "NONE", "BP", 40000)
#mzd_2 = hic.getMatrixZoomData()
print(mzd)
numpy_matrix = mzd.getRecordsAsMatrix(0, 15578282, 0, 15578282)
#numpy_matrix_2 = mzd_2.getRecordsAsMatrix()
plt.imshow(numpy_matrix)
#plt.savefig("/home/ksslab/Siddharth/Temp_plot_dump/HiC_straw_chr_1_0_15178282.png")
plt.show()
plt.close()
#plt.close()
#im = plt.imshow(clr.matrix(balance=False)[:], vmax=2500);
#plt.colorbar(im)
#plt.show()
#plt.close()
def matrix_reader(filename,resolution = 100000,chr_fetch='1',region_fetch =(0,59578282),delete = False):
    if filename[-3:] == "npy":
        return np.load(filename)
    elif filename[-4:] == "cool":
        clr = cooler.Cooler(filename)
        np_mat = clr.matrix(balance=False).fetch((chr_fetch,region_fetch[0],region_fetch[1]))  #59_578_282))
        return np_mat
    elif filename[-3:] == "hic":
        hic = hicstraw.HiCFile(filename)
        mzd = hic.getMatrixZoomData(chr_fetch,chr_fetch,"observed","NONE","BP", resolution)
        numpy_matrix = mzd.getRecordsAsMatrix(region_fetch[0],region_fetch[1],region_fetch[0],region_fetch[1])
        return numpy_matrix
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
    

random_regions_to_check_list = list(zip((str(i) for i in range(1,26,1)),itertools.repeat(15000000), itertools.repeat(30000000)))
#print(random_regions_to_check_list)

#print(os.listdir("/home/ksslab/Siddharth/Program_Directory/Data/Actual_Data/Raw_cool/100kb"))
working_list = list(itertools.product(os.listdir("/home/ksslab/Siddharth/Program_Directory/Data/Actual_Data/Raw_cool/40kb"),random_regions_to_check_list))
print(working_list.index(('2_dpi_S208.cool', ('22', 15000000, 30000000))))
for item in working_list[146:147]:
    print(matrix_reader(f"/home/ksslab/Siddharth/Program_Directory/Data/Actual_Data/Raw_cool/40kb/{item[0]}",40000,item[1][0],(item[1][1],item[1][2])))
    plt.imshow(matrix_reader(f"/home/ksslab/Siddharth/Program_Directory/Data/Actual_Data/Raw_cool/40kb/{item[0]}",40000,item[1][0],(item[1][1],item[1][2])))
    plt.show()
#plt.savefig("/home/ksslab/Siddharth/Temp_plot_dump/Read_mat_chr_1_0_15178282.png")
    #feather.write_feather(pd.DataFrame(matrix_reader(f"/home/ksslab/Siddharth/Program_Directory/Data/Actual_Data/Raw_cool/100kb/{item[0]}",40000,item[1][0],(item[1][1],item[1][2]))),f"/home/ksslab/Siddharth/Program_Directory/Scripts/Python_scripts/Newer/test_runs/100kb/{item[0]}_{item[1][0]}_{item[1][1]}_{item[1][2]}.feather")
#feather.write_feather()
clr2 = cooler.Cooler("/home/ksslab/Siddharth/Temp_Dir/2_dpi_S166.cool")
windows = [2*40000,5*40000]
insulation_table = cooltools.insulation(clr2, windows, verbose=True)
#print(insulation_table)#[['log2_insulation_score_200000']])
#print(feather.read_feather("test.feather").to_numpy()) #Feather is read in as a pandas dataframe like expected so the .to_nnumpy works
#print(np_mat == numpy_matrix,'\n' )
#print(numpy_matrix  == feather.read_feather("test.feather").to_numpy(),'\n')

#plt.imshow(matrix_reader("/home/ksslab/Siddharth/Temp_Dir/2_dpi_S166.cool"))
#plt.show()
#plt.close()

#plt.imshow(matrix_reader("/home/ksslab/Siddharth/Program_Directory/Data/Actual_Data/Raw_cool/40kb/2_dpi_S166.cool"))
#plt.show()
#plt.close()

#print(matrix_reader("/home/ksslab/Siddharth/Temp_Dir/2_dpi_S166.cool") == matrix_reader("/home/ksslab/Siddharth/Program_Directory/Data/Actual_Data/Raw_cool/40kb/2_dpi_S166.cool"))

from cooler.create import ArrayLoader
#import h5py
import cooler
#h = h5py.File("cworld-test_hg19_C-40000-raw.hdf5", 'r')
#heatmap = h['interactions']

# create some bins , using cooler-binnify or some other way
binsize = 40000
chromsizes = pd.read_csv(                                                            
    '/home/ksslab/Siddharth/Program_Directory/Data/Actual_Data/Raw_bams/danRer11.chrom.sizes', 
    sep='\t', 
    names=['name', 'length']).set_index('name')['length']
#print(chromsizes,type(chromsizes),chromsizes.index)
filtered_chromsizes = chromsizes.filter(regex = r'chr1$')#(regex = r'^([a-z]{3}(?:\d{1,2})|X|Y)$')

bins = cooler.binnify(filtered_chromsizes, binsize)
#print(bins)
#print(numpy_matrix.shape)
#help(ArrayLoader)
#print(complete_cooler_mat.shape)
# turn h5oy dataset (2D matrix) into a stream of sparse matrix chunks :
iterator = ArrayLoader(bins, np_mat, chunksize=int(1e6))
#print(iterator)
# load that into cooler:
#help(cooler.create_cooler)
#cooler.create_cooler('output.40kb.cool', bins, iterator, dtypes={"count":"int"}, assembly="dr11")

clr3 = cooler.Cooler('output.40kb.cool')
windows = [2*40000,5*40000]
insulation_table = cooltools.insulation(clr3, windows, verbose=True)
#print(insulation_table)#[['log2_insulation_score_200000']])
#clr4 = cooler.Cooler('yan_temp.cool') #The say S166 hic file converted to cool
clr4 = cooler.Cooler('/home/ksslab/Siddharth/Temp_Dir/again_2_dpi_S208.cool') #The say S208 hic file converted to cool
windows = [2*40000, 5*40000]
insulation_table = cooltools.insulation(clr4,windows, verbose=True)
print(clr4.matrix)
np_mat = clr4.matrix(balance=False).fetch(('22',15000000,30000000))  #59_578_282))
hmap = plt.imshow(np_mat,cmap = 'Spectral')
plt.colorbar(hmap)
plt.show()
plt.close()
print(np_mat)
#print(insulation_table.columns.tolist())
#print(set(insulation_table['log2_insulation_score_80000']))
#print(insulation_table[insulation_table['log2_insulation_score_200000'].notnull()])
#print(pd.DataFrame(insulation_table['log2_insulation_score_200000']),insulation_table['log2_insulation_score_200000'])
#feather.write_feather(pd.DataFrame(insulation_table[['log2_insulation_score_200000','log2_insulation_score_80000']]),'test_insul_df.feather')
#np.save("size_test_mat_2.npy",np_mat)
clr = cooler.Cooler("/home/ksslab/Siddharth/Program_Directory/Data/Actual_Data/Raw_cool/40kb/2_dpi_S208.cool")
np_mat = clr.matrix(balance=False).fetch(('22',15000000,30000000))  #59_578_282))
hmap = plt.imshow(np_mat,cmap = 'Spectral')
plt.colorbar(hmap)
plt.show()
plt.close()
print(np_mat,'\n\n')
        

clr = cooler.Cooler("/home/ksslab/Siddharth/Temp_Dir/Temp_cooler_S208_2_dpi.cool")
np_mat = clr.matrix(balance=False).fetch(('chr22',15000000,30000000))  #59_578_282))
hmap = plt.imshow(np_mat,cmap = 'Spectral')
plt.colorbar(hmap)
plt.show()
plt.close()
print(np_mat)

dr11_chromsizes = bioframe.fetch_chromsizes("danRer11")
#print(dr11_chromsizes)
#print(bioframe.assembly_info('danRer11'))
#dr11_cens = bioframe.fetch_centromeres('hg38',provider='ucsc')
# create a view with chromosome arms using chromosome sizes and definition of centromeres
#dr11_arms = bioframe.make_chromarms(dr11_chromsizes,  dr11_cens)
#print(dr11_arms)
# select only those chromosomes available in cooler
#dr11_arms_arms = dr11_arms[dr11_arms.chrom.isin(clr.chromnames)].reset_index(drop=True)
#print(dr11_arms)