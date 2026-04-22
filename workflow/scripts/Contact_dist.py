# import core packages
import warnings
warnings.filterwarnings("ignore")
from itertools import combinations
import os

# import semi-core packages
import matplotlib.pyplot as plt
from matplotlib import colors

plt.style.use('seaborn-v0_8-poster')
import numpy as np
import pandas as pd
from multiprocessing import Pool

# import open2c libraries
import bioframe

import cooler
import cooltools

from packaging import version
if version.parse(cooltools.__version__) < version.parse('0.5.2'):
    raise AssertionError("tutorial relies on cooltools version 0.5.2 or higher,"+
                         "please check your cooltools version and update to the latest")

# count cpus
num_cpus = os.getenv('SLURM_CPUS_PER_TASK')
if not num_cpus:
    num_cpus = os.cpu_count()
num_cpus = int(num_cpus)



# Load a Hi-C map at a 1kb resolution from a cooler file.
resolution = 50000 # note this might be slightly slow on a laptop
                  # and could be lowered to 10kb for increased speed
possibilities = ['control_S68','12_hpi_S207','2_dpi_S208','4_dpi_S69','7_dpi_S1']
name_color_dict= {'control_S68':'purple','12_hpi_S207':'blue','2_dpi_S208':'green','4_dpi_S69':'orangered','7_dpi_S1':'red'}
f, ax = plt.subplots(1,1)
for name in possibilities:
    clr = cooler.Cooler(f'/home/ksslab/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/temporaries/cools/{name}.mcool::/resolutions/{resolution}')
    cvd_smooth_agg = cooltools.expected_cis(
        clr=clr,
        #view_df=hg38_arms,
        smooth=True,
        aggregate_smoothed=True,
        smooth_sigma=0.1,
        nproc=num_cpus
    )
    cvd_smooth_agg['balanced.avg.smoothed.agg'].loc[cvd_smooth_agg['dist'] < 2] = np.nan

    for region in clr.chromnames:
        if region == clr.chromnames[len(clr.chromnames)-1]:
            ax.loglog(
            cvd_smooth_agg['dist_bp'].loc[cvd_smooth_agg['region1']==region],
            cvd_smooth_agg['balanced.avg.smoothed.agg'].loc[cvd_smooth_agg['region1']==region],label = f"{name}",color = name_color_dict[name],linewidth =1 
        )
        else:
            ax.loglog(
                cvd_smooth_agg['dist_bp'].loc[cvd_smooth_agg['region1']==region],
                cvd_smooth_agg['balanced.avg.smoothed.agg'].loc[cvd_smooth_agg['region1']==region],color = name_color_dict[name],linewidth =1 
             )
        ax.set(
            xlabel='separation, bp',
            ylabel='IC contact frequency')
        ax.set_aspect(1.0)
        ax.grid(lw=0.5)
plt.title(f"Contact Frequency vs Distance for various timepoints at {resolution} bin size")
plt.legend()
plt.savefig(f"/home/ksslab/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/contact_freq/{resolution}/contact_freq_dist_cis_{resolution}.svg")
plt.show()

'''
cvd = cooltools.expected_cis(
    clr=clr,
    #view_df=hg38_arms,
    smooth=False,
    aggregate_smoothed=False,
    nproc=num_cpus #if you do not have multiple cores available, set to 1
)
print(cvd)

f, ax = plt.subplots(1,1)

for region in clr.chromnames:
    ax.loglog(
        cvd['dist_bp'].loc[cvd['region1']==region],
        cvd['contact_frequency'].loc[cvd['region1']==region],
    )
    ax.set(
        xlabel='separation, bp',
        ylabel='IC contact frequency')
    ax.set_aspect(1.0)
    ax.grid(lw=0.5)
plt.show()
'''
    

#cvd_smooth_agg['balanced.avg.smoothed'].loc[cvd_smooth_agg['dist'] < 2] = np.nan

#f, ax = plt.subplots(1,1)

#for region in clr.chromnames:
    #ax.loglog(
        #cvd_smooth_agg['dist_bp'].loc[cvd_smooth_agg['region1']==region],
        #cvd_smooth_agg['balanced.avg.smoothed'].loc[cvd_smooth_agg['region1']==region],
    #)
    #ax.set(
        #xlabel='separation, bp',
        #ylabel='IC contact frequency')
    #ax.set_aspect(1.0)
    #ax.grid(lw=0.5)
#plt.show()


