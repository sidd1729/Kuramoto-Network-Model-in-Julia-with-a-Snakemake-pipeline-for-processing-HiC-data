
import cooler
import cooltools
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyarrow.feather as feather
import os
import h5py
import bioframe
import sys
from sklearn.decomposition import PCA #Performing principal component analysis of HiC matrices
from tqdm import tqdm

#System arguments are 1: Cool file path name, 2: Resolution, 3: Window size
clr = cooler.Cooler(f"{sys.argv[1]}::/resolutions/{sys.argv[2]}")
#clr2 = cooler.balance_cooler(clr)
window_argument = int(sys.argv[3])*int(sys.argv[2])
windows = [window_argument]

#Insulation calculation
insulation_table = cooltools.insulation(clr, windows, verbose=True)

#Expected value calculation
num_cpus =1
cvd = cooltools.expected_cis(
    clr=clr,
    #observed=False,
    #view_df=hg38_arms,
    smooth=True,
    aggregate_smoothed=True,
    smooth_sigma=0.1,
    nproc=num_cpus #if you do not have multiple cores available, set to 1
)

#print(insulation_table,insulation_table[f'log2_insulation_score_{window_argument}'])
#print(cvd.columns)
#print(cvd[[ 'count.sum', 'balanced.sum', 'count.avg', 'balanced.avg']])
#print(cvd)

feather.write_feather(pd.DataFrame(insulation_table[f'log2_insulation_score_{window_argument}']),f"{sys.argv[4]}/{os.path.splitext(os.path.basename(sys.argv[1]))[0]}_{sys.argv[2]}_insul.feather")

for i in tqdm(clr.chromnames):
    np_mat = clr.matrix(balance=False).fetch(f'{i}')
    insulation_table_chrom = insulation_table[(insulation_table['chrom'] == i) & (insulation_table['region'] == i)]
    #print(f"{insulation_table_chrom} is insulation table for specific chromosome")
    #print(f"{insulation_table_chrom[f'log2_insulation_score_{window_argument}']} is insulation score for specific chromosome")
    feather.write_feather(pd.DataFrame(insulation_table_chrom[f'log2_insulation_score_{window_argument}']),f"{sys.argv[4]}/{os.path.splitext(os.path.basename(sys.argv[1]))[0]}_{sys.argv[2]}_{i}_insul.feather")
    #hmap = plt.imshow(np_mat,cmap = 'Spectral_r')
    #plt.colorbar(hmap)
    #plt.show()
    #plt.close()

    #Generating an expected matrix
    trial = np.zeros(np_mat.shape)
    cvd_tl = cvd[(cvd['region1'] == f'{i}') & (cvd['region2'] == f'{i}')]

    rows, cols = np.indices(trial.shape)
    diffs = np.abs(rows - cols)
    for dist in cvd_tl['dist']:
        mask = diffs == dist  # positions where abs(i-j) == dist
        #print(mask)
        cvd_tl_work = cvd_tl[cvd_tl['dist'] == dist]
        #print(cvd_tl_work['count.avg'])
        count_avg_len = cvd_tl_work['count.avg']
        #print(len(count_avg_len),count_avg_len)
        #print(trial[mask])
        trial[mask] = count_avg_len

    feather.write_feather(pd.DataFrame(trial),f"{sys.argv[5]}/{os.path.splitext(os.path.basename(sys.argv[1]))[0]}_{sys.argv[2]}_{i}_exp.feather")
    #hmap = plt.imshow(trial,cmap = 'Spectral_r')
    #plt.colorbar(hmap)
    #plt.title("Raw Expected Chr1 2dpi S166 50kb")
    #plt.savefig("Raw_Exp_1_2dpi_S166_50kb.png")
    #plt.show()
    #plt.close()
    
    o_min_e = np_mat-trial
    feather.write_feather(pd.DataFrame(o_min_e),f"{sys.argv[6]}/{os.path.splitext(os.path.basename(sys.argv[1]))[0]}_{sys.argv[2]}_{i}_o_min_e.feather")
    #hmap = plt.imshow(np.log(o_min_e),cmap = 'Spectral_r')
    #plt.colorbar(hmap)
    #plt.title("Log Observed-Expected Chr1 2dpi S166 50kb")
    #plt.savefig("Log_O_min_E_1_2dpi_S166_50kb.png")
    #plt.show()
    #plt.close()

    o_by_e = np.divide(np_mat,trial)
    feather.write_feather(pd.DataFrame(o_by_e),f"{sys.argv[7]}/{os.path.splitext(os.path.basename(sys.argv[1]))[0]}_{sys.argv[2]}_{i}_o_by_e.feather")
    #hmap = plt.imshow(np.log(o_by_e),cmap = 'Spectral_r')
    #plt.colorbar(hmap)
    #plt.show()
    #plt.close()
    
    '''
    #PCs
    no_components = 3
    pca = PCA(n_components = no_components)
    if not np.any(np.isnan(o_min_e)):
        PCA_output = pca.fit(o_min_e)
        ICEd_PCs = PCA_output.components_
        plt.plot(list(range(len(ICEd_PCs[0,:]))),ICEd_PCs[0,:])
        plt.show()
        plt.close()
        plt.plot(list(range(len(ICEd_PCs[1,:]))),ICEd_PCs[1,:])
        plt.show()
        plt.close()
        plt.plot(list(range(len(ICEd_PCs[2,:]))),ICEd_PCs[2,:])
        plt.show()
        plt.close()
    
    if not np.any(np.isnan(o_by_e)):
        PCA_output = pca.fit(o_by_e)
        ICEd_PCs = PCA_output.components_
        plt.plot(list(range(len(ICEd_PCs[0,:]))),ICEd_PCs[0,:])
        plt.show()
        plt.close()
        plt.plot(list(range(len(ICEd_PCs[1,:]))),ICEd_PCs[1,:])
        plt.show()
        plt.close()
        plt.plot(list(range(len(ICEd_PCs[2,:]))),ICEd_PCs[2,:])
        plt.show()
        plt.close()
    '''

