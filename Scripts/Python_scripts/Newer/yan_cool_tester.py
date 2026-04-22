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
import fanc
import time
import sys
#print(sys.argv[0],sys.argv[1])
clr = cooler.Cooler("/home/ksslab/Siddharth/Program_Directory/Data/Actual_Data/Raw_cool/40kb/control_S63.cool")
chromstarts = []
for i in clr.chromnames:
    print(f'{i} : {clr.extent(i)}')
    chromstarts.append(clr.extent(i)[0])
np_mat = clr.matrix(balance=False).fetch('1')
print(np_mat,'\n\n',np.sum(np_mat))
#feather.write_feather(pd.DataFrame(np_mat),"/home/ksslab/Siddharth/Program_Directory/Data/Intermediate_Processing_Data/control_S63.feather")

clr2 = cooler.Cooler("/home/ksslab/Siddharth/Program_Directory/Data/Actual_Data/Raw_cool/40kb/12_hpi_S165.cool")
chromstarts = []
for i in clr2.chromnames:
    print(f'{i} : {clr2.extent(i)}')
    chromstarts.append(clr2.extent(i)[0])
np_mat = clr2.matrix(balance=False).fetch('1')
print(np_mat,'\n\n',np.sum(np_mat))
#feather.write_feather(pd.DataFrame(np_mat),"/home/ksslab/Siddharth/Program_Directory/Data/Intermediate_Processing_Data/12_hpi_S165.feather")

clr3 = cooler.Cooler("/home/ksslab/Siddharth/2dpi.mcool::resolutions/50000")
chromstarts = []
for i in clr3.chromnames:
    print(f'{i} : {clr3.extent(i)}')
    chromstarts.append(clr3.extent(i)[0])
np_mat = clr3.matrix(balance=False).fetch(('22',15000000,30000000))  #59_578_282))
hmap = plt.imshow(np_mat,cmap = 'Spectral_r')
plt.colorbar(hmap)
plt.show()
plt.close()
print(np_mat,'\n\n')

#hic_500kb = fanc.load("/home/ksslab/Siddharth/Program_Directory/Data/Actual_Data/Raw_cool/40kb/2_dpi_S208.cool")
#intra_expected, intra_expected_chromosome, inter_expected = hic_500kb.expected_values()
#print(intra_expected, intra_expected_chromosome, inter_expected)


# obtain bin distances
#bin_size = hic_500kb.bin_size
#distance = list(range(0, bin_size * len(intra_expected_chromosome['1']), bin_size))

# plot expected values
#fig, ax = plt.subplots()
#plt.plot(distance, intra_expected_chromosome['1'])
#ax.set_xscale('log')
#ax.set_yscale('log')
#ax.set_xlabel("Distance")
#ax.set_ylabel("Average contacts")
#plt.show()

num_cpus = 1
start_time = time.perf_counter()
'''
cvd = cooltools.expected_cis(
    clr=clr,
    #view_df=hg38_arms,
    smooth=True,
    aggregate_smoothed=True,
    smooth_sigma=0.1,
    nproc=num_cpus #if you do not have multiple cores available, set to 1
)
'''
end_time = time.perf_counter()
print(f"{end_time-start_time} is time elapsed")
'''
print(cvd.head(4))
print(cvd.tail(4))
print(type(cvd))
'''