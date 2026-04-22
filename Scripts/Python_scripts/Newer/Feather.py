import pyarrow.feather as feather
import pandas as pd
import pickle
def general_unpickler(filename):
    with open(filename,'rb') as f:
        unp_object = pickle.load(f)
        return unp_object
    
GM12878_whole_chr1_kura_r = general_unpickler("/home/ksslab/Siddharth/Program_Directory/Data/Simulation_Results/Linux/GM12878/Chr1_Kura_r.pickle")
feather.write_feather(pd.DataFrame(GM12878_whole_chr1_kura_r ),"/home/ksslab/Siddharth/Program_Directory/Data/Simulation_Results/Linux/GM12878/Chr1_Kura_r.feather")