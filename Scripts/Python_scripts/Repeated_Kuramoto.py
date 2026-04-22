import numpy as np
import csv
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