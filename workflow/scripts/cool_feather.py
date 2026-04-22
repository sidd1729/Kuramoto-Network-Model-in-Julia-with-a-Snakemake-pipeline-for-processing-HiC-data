import cooler
import cooltools
import matplotlib.pyplot as plt
import hicstraw
import numpy as np
import pandas as pd
import pyarrow.feather as feather
import os
from tqdm import tqdm
import sys
#print(sys.argv[0],sys.argv[1])
clr = cooler.Cooler(f"{sys.argv[1]}::resolutions/{sys.argv[2]}")
#clr = cooler.Cooler("/home/ksslab/Siddharth/Program_Directory/Data/Actual_Data/Raw_cool/40kb/control_S63.cool")
#chromstarts = []
for i in tqdm(clr.chromnames):
    #print(f'{i} : {clr.extent(i)}')
    #chromstarts.append(clr.extent(i)[0])
    np_mat = clr.matrix(balance=False).fetch(i)
    #print(np_mat,'\n\n',np.sum(np_mat),i)
    feather.write_feather(pd.DataFrame(np_mat),f"{sys.argv[3]}/{os.path.splitext(os.path.basename(sys.argv[1]))[0]}_{sys.argv[2]}_{i}.feather")
 
