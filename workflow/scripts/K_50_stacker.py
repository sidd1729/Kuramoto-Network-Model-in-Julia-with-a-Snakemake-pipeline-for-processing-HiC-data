import pyarrow.feather as feather
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
import cmasher
import os
column_list = []
sample_names = ["control_S68","12_hpi_S207","2_dpi_S208","4_dpi_S69","7_dpi_S1"]
chromosome_list = [str(x) for x in range(1,26)]
dir_sorted  = os.listdir("/home/ksslab/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/50000/k_50")
filtered_control = sorted(filter(lambda x: (f"{x.split('_')[0]}_{x.split('_')[1]}" == "control_S68") & (f"{x.split('_')[4]}" == "kura"),dir_sorted),key = lambda x : int(f"{x.split('_')[3]}"))
filtered_12_hpi = sorted(filter(lambda x: (f"{x.split('_')[0]}_{x.split('_')[1]}" == "12_hpi") & (f"{x.split('_')[5]}" == "kura"),dir_sorted),key = lambda x : int(f"{x.split('_')[4]}"))
filtered_2_dpi = sorted(filter(lambda x: (f"{x.split('_')[0]}_{x.split('_')[1]}" == "2_dpi") & (f"{x.split('_')[5]}" == "kura"),dir_sorted),key = lambda x : int(f"{x.split('_')[4]}"))
filtered_4_dpi = sorted(filter(lambda x: (f"{x.split('_')[0]}_{x.split('_')[1]}" == "4_dpi") & (f"{x.split('_')[5]}" == "kura"),dir_sorted),key = lambda x : int(f"{x.split('_')[4]}"))
filtered_7_dpi = sorted(filter(lambda x: (f"{x.split('_')[0]}_{x.split('_')[1]}" == "7_dpi") & (f"{x.split('_')[5]}" == "kura"),dir_sorted),key = lambda x : int(f"{x.split('_')[4]}"))
#print(f"{list(filtered_control)} \n {list(filtered_12_hpi)} \n {list(filtered_2_dpi)} \n {list(filtered_4_dpi)} \n {list(filtered_7_dpi)}")

chromosome_bin_size  = []
for file in filtered_control:
    chromosome_bin_size.append(feather.read_feather(f"/home/ksslab/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/50000/k_50/{file}").shape[0])
#print(chromosome_bin_size)
chromosome_size_dict = {'bins for chromosomes' : chromosome_bin_size}
chromosome_size_df = pd.DataFrame(chromosome_size_dict)
feather.write_feather(chromosome_size_df,"/home/ksslab/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/k_50_table/k_50_chrom_bin_nos.feather")
#map(lambda x: feather.read_feather(f"/home/ksslab/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/50000/k_50/{x}).to_numpy()")

def k_50_value_extractor(k_50_df):
    k_50_values = map(lambda x : x[0][0][0],np.array(k_50_df))
    return k_50_values


#print((feather.read_feather(f"/home/ksslab/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/50000/k_50/control_S68_50000_1_kura_k_50.feather").to_numpy()))
#print("\n\n")
control_list  = []
for control_file in filtered_control:
    control_list.extend(list(k_50_value_extractor(feather.read_feather(f"/home/ksslab/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/50000/k_50/{control_file}"))))
#print(len(control_list),"\n",control_list[0:4])
hpi_12_list  = []
for hpi_12_file in filtered_12_hpi:
    hpi_12_list.extend(list(k_50_value_extractor(feather.read_feather(f"/home/ksslab/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/50000/k_50/{hpi_12_file}"))))
#print(len(hpi_12_list),"\n",hpi_12_list[0:4])
dpi_2_list  = []
for dpi_2_file in filtered_2_dpi:
    dpi_2_list.extend(list(k_50_value_extractor(feather.read_feather(f"/home/ksslab/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/50000/k_50/{dpi_2_file}"))))
#print(len(dpi_2_list),"\n",dpi_2_list[0:4])
dpi_4_list  = []
for dpi_4_file in filtered_4_dpi:
    dpi_4_list.extend(list(k_50_value_extractor(feather.read_feather(f"/home/ksslab/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/50000/k_50/{dpi_4_file}"))))
#print(len(dpi_4_list),"\n",dpi_4_list[0:4])
dpi_7_list  = []
for dpi_7_file in filtered_7_dpi:
    dpi_7_list.extend(list(k_50_value_extractor(feather.read_feather(f"/home/ksslab/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/50000/k_50/{dpi_7_file}"))))
#print(len(dpi_7_list),"\n",dpi_7_list[0:4])

data_dict = {"control" : control_list, "12 hpi" : hpi_12_list, "2 dpi" : dpi_2_list , "4 dpi" : dpi_4_list, "7 dpi" : dpi_7_list }
#print(pd.DataFrame(data_dict),pd.DataFrame(data_dict).to_numpy())
data_frame = pd.DataFrame(data_dict)
feather.write_feather(data_frame,"/home/ksslab/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/k_50_table/k_50_sample_wise_df.feather")
sample_labels = ["control","12 hpi","2 dpi ","4 dpi","7 dpi"]
np_plot = data_frame.to_numpy()
#print(np.max(np_plot))
cmap = plt.get_cmap("viridis")
cmap.set_under("black")
plt.figure()
#sns.heatmap(np_plot)
plt.imshow(np_plot, aspect="auto",cmap = cmap,vmin=0.0,vmax=np.max(np_plot))
plt.title("Heatmap of K_50 values samplewise")
plt.colorbar()
plt.xticks(ticks = np.arange(5),labels = sample_labels,rotation = 45)
plt.savefig("/home/ksslab/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/k_50_table/k_50_heatmap.svg")
plt.show()
#help(sorted)
#print((lambda x : f"{x.split('_')[4]}")("4_dpi_S69_50000_2_kura_k_50.feather"))

#Sort using the chromosome number as key from the split
#print(feather.read_feather("/home/ksslab/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/k_50_table/k_50_chrom_bin_nos.feather"))
print(feather.read_feather("/home/ksslab/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/k_50_table/k_50_sample_wise_df.feather"))