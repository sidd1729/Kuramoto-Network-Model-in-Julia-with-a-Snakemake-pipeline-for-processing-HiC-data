import pyarrow.feather as feather
import numpy as np
import matplotlib.pyplot as plt
import os
from scipy.stats import norm,gaussian_kde,kde
from sklearn.neighbors import KernelDensity
from tqdm import tqdm

#chromosome_list = [str(x) for x in range(1,26)]
#empty_list = []
#trial_attempt = list(map(lambda x: x.split('_')[3:5], os.listdir("/home/ksslab/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/PCs/50000")))
#print(trial_attempt)
#bins= np.arange(-3,4,1) 
range_to_plot = np.arange(-5,5)
plt.style.use('bmh')
resolution = 10000

#file = "4DNFITRVKRPA_50000_PC1.feather"
#for i in chromosome_list:
plt.figure()
for file in sorted(os.listdir(f"/home/ksslab/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/PCs/{resolution}")):
    #if ((file.split('_')[4] == i) | (file.split('_')[3] == i)) & (file[-11:-8:+1] == "PC1"):
    if (file[-11:-8:+1] == "PC1") & ((file.split('_'))[-2] == f'{resolution}'):
        print(file)
#file = f"12_hpi_S207_{resolution}_PC1.feather"
        trial = np.array(feather.read_feather(f"/home/ksslab/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/PCs/{resolution}/{file}"))
#trial = np.array(feather.read_feather("/home/ksslab/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/temporaries/trials/4DNFITRVKRPA_50000_PC1.feather"))
        #counts, bins = np.histogram(trial)
#print(trial,trial.dtypes)
#print(np.isfinite(trial))
#finite_trial = trial.dropna().to_numpy()
        finite_trial = trial[np.isfinite(trial)] #Extracting only finite elements from input array
#print(finite_trial)
#kde = KernelDensity(kernel='gaussian', bandwidth=0.5).fit(finite_trial.reshape(1,-1))
#density =kde.score_samples(finite_trial)
        test = np.linspace(-5,5,len(finite_trial))
        density = gaussian_kde(finite_trial)
#print(density(test))
        plt.plot(test,density(test),label = f"{file.split('_')[0]}_{file.split('_')[1]}",linewidth = 1)
        #normed = list(map(lambda x: x/np.nansum(np.array(trial)),np.array(trial)))
        #plt.plot(list(range(len(normed))),normed, label = f"{file.split('_')[0]}_{file.split('_')[1]}")
        #plt.hist(np.array(trial),bins = 20,density=True,alpha =0.8,label = f"{file.split('_')[0]}_{file.split('_')[1]}")
plt.title(f"PC1 probability density for for {int(resolution/1000)} kb bins")
plt.xlabel("PC1 value")
plt.ylabel("Probability density")
plt.legend()
plt.savefig(f"/home/ksslab/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/Dists/PC1/{resolution}/PC1_dist.png")
plt.show()

    #plt.figure()
    #for file in sorted(os.listdir(f"/home/ksslab/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/insulation/{resolution}")):
        #if ((file.split('_')[4] == i) | (file.split('_')[3] == i)):
            #print(file)
            #trial_2 = feather.read_feather(f"/home/ksslab/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/insulation/{resolution}/{file}")
            #normed = list(map(lambda x: x/np.nansum(np.array(trial_2)),np.array(trial_2)))
            #plt.plot(list(range(len(normed))),normed, label = f"{file.split('_')[0]}_{file.split('_')[1]}")
            #plt.hist(np.array(trial_2),bins = 20,density=True,alpha =0.8,label = f"{file.split('_')[0]}_{file.split('_')[1]}")
    #plt.title(f"Insulation probability density for chr{i} for {int(resolution/1000)} kb bins")
    #plt.xlabel("Insulation score")
    #plt.ylabel("Probability density")
    #plt.legend()
    #plt.savefig(f"/home/ksslab/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/Dists/insul/{resolution}/insul_dist_chr_{i}.png")
    #plt.show()

#trial_2 = feather.read_feather("/home/ksslab/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/PCs/50000/2_dpi_S208_50000_1_PC1.feather")
#normed_2 = list(map(lambda x: x/np.nansum(trial_2.iloc[:,0]),np.array(trial_2)))

#plt.plot(list(range(len(normed_2))),normed_2)
#plt.show()
