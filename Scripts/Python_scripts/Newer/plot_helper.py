import numpy as np
import pandas as pd
import csv
import pyarrow.feather as feather
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
    

feather.write_feather(pd.DataFrame(matrix_reader("/home/ksslab/Siddharth/Program_Directory/Data/hg38/GM12878/100kb/raw/chr1.mat.txt")),"/home/ksslab/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/temporaries/trials/GM128278_100kb_chr1.feather")