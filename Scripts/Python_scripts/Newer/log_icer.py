import cmasher
from iced import normalization
import numpy as np
import os
import matplotlib.pyplot as plt
import pyarrow.feather as feather

for item in os.listdir("/home/ksslab/Siddharth/Program_Directory/Scripts/Python_scripts/Newer/test_runs/100kb"):
    print(item)
    relevant_mat = np.log(normalization.ICE_normalization(feather.read_feather(f"/home/ksslab/Siddharth/Program_Directory/Scripts/Python_scripts/Newer/test_runs//100kb/{item}").to_numpy()) +1)
    hmap = plt.imshow(relevant_mat,cmap = 'cmr.sunburst_r')
    plt.colorbar(hmap)
    plt.title(f"{item}")
    plt.savefig(f"/home/ksslab/Siddharth/Temp_plot_dump/100kb_Real/{item}.png")
    #plt.show()
    plt.close()