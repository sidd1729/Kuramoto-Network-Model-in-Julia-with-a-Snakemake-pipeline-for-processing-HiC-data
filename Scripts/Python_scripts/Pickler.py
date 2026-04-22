import pickle
import os
import sys
def general_unpickler(filename):
    with open(filename,'rb') as f:
        unp_object = pickle.load(f)
        return unp_object


def general_pickler(filename,tbp_object=None,read = False,reader_fn = None):
    if read == True:
        tbp_object = reader_fn(filename)      
    with open(filename + ".pickle",'wb') as f:
        pickle.dump(tbp_object,f)

def tuple_list_reader(filename):
    file = open('{}'.format(filename), 'r')
    content = file.read().splitlines()
    tuple_list =[]
    file.close()
    for row in content:
        tuple_list.append(tuple(float(i) for i in row.split(' ')))
    return tuple_list

dir_path = sys.argv[1] #Directory containing files to be pickled
dir_list = os.listdir(dir_path)
print(dir_list)
for file in dir_list:
    if file[-7:] != ".pickle":
        general_pickler(dir_path + '/'+ file,read = True,reader_fn=tuple_list_reader)