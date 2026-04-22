import subprocess
import os
import pickle
#dir_path = "C:/Users/HP/Downloads/Siddharth/Scripts/Python_scripts/HESC"
#dir_list = os.listdir(dir_path)
#print(dir_list)

def general_pickler(filename,tbp_object=None,read = False,reader_fn = None):
    if read == True:
        tbp_object = reader_fn(filename)      
    with open(filename + ".pickle",'wb') as f:
        pickle.dump(tbp_object,f)
'''
file_list = []
for file in dir_list:
    if file[-7:] != ".pickle":
        #subprocess.call(['python', 'Kuramoto_example.py','system',"C:/Users/HP/Downloads/Siddharth/Data/hg38/GM12878/100kb/raw/chr1.mat.txt","(143000001,158000001)",100000,10,"C:/Users/HP/Downloads/Siddharth/Data/hg38/GM12878/100kb/Chromosome_1","C:/Users/HP/Downloads/Siddharth/Thesis_Figures/Weighted_Network/GM12878/100kb/Chromosome_1"])
        file_list.append(file)
print(file_list)
'''
trial_list = []
for i in range(0,146): #0,141 for 10 x 10
    try:
        trial = subprocess.check_output(['python','R_python_interface.py',"/home/ksslab/Siddharth/Program_Directory/Data/Simulation_Results/Windows/HESC_clean/ICE_5_5/10k_coupling_k_vals__500_Chr1_100kb_HESC_ICE_5_5_first_15_{}.pickle".format(i),"/home/ksslab/Siddharth/Program_Directory/Data/Simulation_Results/Windows/HESC_clean/ICE_5_5/10k_coupling_trial_500_Chr1_100kb__HESC_ICE_5_5_first_15_{}.pickle".format(i),"/home/ksslab/Siddharth/Temp_plot_dump/plot_{}_HESC_ICE_5_5.png".format(i),"k-r for Raw-Kura order chr1 100kb index {} of 5x5 HESC ICE".format(i)],text=True)
        #print(type(trial),trial,trial[:len(trial-1)],float(trial[:len(trial-1)]))
        #print(trial , type(trial), trial[0:len(trial)-1]) #tuple(map(float,trial[0:len(trial)-1])))
        #print(trial[0:len(trial)-1].split(',')))
        floated = tuple(map(float,trial[0:len(trial)-1].split(',')))
        print(floated)
        #print(floated,type(floated),type(floated[0]),type(floated[1]))
        trial_list.append((floated,i))
    except:
        continue
print(trial_list)
general_pickler("/home/ksslab/Siddharth/Program_Directory/Data/Simulation_Results/Windows/HESC_ICE_K_50_5_5_trial_list.pickle",trial_list)

#Redo
        