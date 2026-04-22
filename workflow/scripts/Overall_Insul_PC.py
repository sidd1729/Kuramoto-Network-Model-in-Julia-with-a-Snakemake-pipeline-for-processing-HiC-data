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

resolution = 500000 #Default will be 50000 if not specified in output filename
window_size = 5
#System arguments are 1: Cool file path name, 2: Resolution, 3: Window size
clr1 = cooler.Cooler(f"/home/ksslab/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/temporaries/cools/control_S68.mcool::resolutions/{resolution}")
clr3 = cooler.Cooler(f"/home/ksslab/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/temporaries/cools/2_dpi_S208.mcool::resolutions/{resolution}")
clr2 = cooler.Cooler(f"/home/ksslab/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/temporaries/cools/12_hpi_S207.mcool::resolutions/{resolution}")
clr4 = cooler.Cooler(f"/home/ksslab/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/temporaries/cools/4_dpi_S69.mcool::resolutions/{resolution}")
clr5 = cooler.Cooler(f"/home/ksslab/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/temporaries/cools/7_dpi_S1.mcool::resolutions/{resolution}")

#clr2 = cooler.balance_cooler(clr)
#window_argument = int(sys.argv[3])*int(sys.argv[2])
window_argument_for_now = int(window_size*resolution)
windows = [window_argument_for_now]

#Insulation calculation
insulation_table_1 = cooltools.insulation(clr1, windows, verbose=True)
insulation_table_2 = cooltools.insulation(clr2, windows, verbose=True)
insulation_table_3 = cooltools.insulation(clr3, windows, verbose=True)
insulation_table_4 = cooltools.insulation(clr4, windows, verbose=True)
insulation_table_5 = cooltools.insulation(clr5, windows, verbose=True)
insulation_timepoint_df = pd.concat([insulation_table_1['chrom'],insulation_table_1['start'],insulation_table_1['end'],insulation_table_1[f'log2_insulation_score_{window_argument_for_now}'],insulation_table_2[f'log2_insulation_score_{window_argument_for_now}'],insulation_table_3[f'log2_insulation_score_{window_argument_for_now}'],insulation_table_4[f'log2_insulation_score_{window_argument_for_now}'],insulation_table_5[f'log2_insulation_score_{window_argument_for_now}']],axis =1 ,keys = ['chrom','start','end','control','12hpi','2dpi','4dpi','7dpi'])
insulation_timepoint_df['Chomosome_Region'] = insulation_timepoint_df.apply(lambda row : f"chr{row['chrom']}:{row['start']}-{row['end']}", axis = 1)
#insulation_timepoint_df_pos = insulation_timepoint_df.pop('pos')

print(insulation_timepoint_df)
#feather.write_feather(pd.DataFrame(insulation_table[f'log2_insulation_score_{window_argument}']),f"{sys.argv[4]}/{os.path.splitext(os.path.basename(sys.argv[1]))[0]}_{sys.argv[2]}_{i}_insul.feather")
feather.write_feather(insulation_timepoint_df,f"/home/ksslab/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/Overall_table/Overall_insulation_{resolution}.feather")

bins_1 = clr1.bins()[:]
bins_2 = clr2.bins()[:]
bins_3 = clr3.bins()[:]
bins_4 = clr4.bins()[:]
bins_5 = clr5.bins()[:]

dr_11_genome = bioframe.load_fasta("/home/ksslab/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/temporaries/genome/danRer11.fa")

#print(dr_11_genome.items())
dr_11_genome_renamed = {chrom.replace('chr', '') if chrom.startswith('chr') else chrom: seq for chrom, seq in dr_11_genome.items()}
#print(dr_11_genome.items())
#print(bins)
#print(type(bins['chrom']),type(map(lambda x: f'chr{x}',bins['chrom'])))
#def chrom_converter(x):
    #return f'chr{x}'
#bins['chrom'] = bins['chrom'].map(chrom_converter)
#print(bins)
gc_cov_1 = bioframe.frac_gc(bins_1[['chrom','start','end']],dr_11_genome_renamed) #genome renamed to be used
gc_cov_2 = bioframe.frac_gc(bins_2[['chrom','start','end']],dr_11_genome_renamed) #genome renamed to be used
gc_cov_3 = bioframe.frac_gc(bins_3[['chrom','start','end']],dr_11_genome_renamed) #genome renamed to be used
gc_cov_4 = bioframe.frac_gc(bins_4[['chrom','start','end']],dr_11_genome_renamed) #genome renamed to be used
gc_cov_5 = bioframe.frac_gc(bins_5[['chrom','start','end']],dr_11_genome_renamed) #genome renamed to be used
#print(gc_cov)
#chrom_name_values = list(map(chrom_converter,clr.chromnames))
#print(type(clr.chromnames))
view_df_1 = pd.DataFrame({'chrom': clr1.chromnames,
                        'start': 0,
                        'end': clr1.chromsizes.values,
                        'name': clr1.chromnames})

view_df_2 = pd.DataFrame({'chrom': clr2.chromnames,
                        'start': 0,
                        'end': clr2.chromsizes.values,
                        'name': clr2.chromnames})

view_df_3 = pd.DataFrame({'chrom': clr3.chromnames,
                        'start': 0,
                        'end': clr3.chromsizes.values,
                        'name': clr3.chromnames})

view_df_4 = pd.DataFrame({'chrom': clr4.chromnames,
                        'start': 0,
                        'end': clr4.chromsizes.values,
                        'name': clr4.chromnames})

view_df_5 = pd.DataFrame({'chrom': clr5.chromnames,
                        'start': 0,
                        'end': clr5.chromsizes.values,
                        'name': clr5.chromnames})

#print(view_df)
cis_eigs_1 = cooltools.eigs_cis(clr1,gc_cov_1,view_df=view_df_1,n_eigs=3)
cis_eigs_2 = cooltools.eigs_cis(clr2,gc_cov_2,view_df=view_df_2,n_eigs=3)
cis_eigs_3 = cooltools.eigs_cis(clr3,gc_cov_3,view_df=view_df_3,n_eigs=3)
cis_eigs_4 = cooltools.eigs_cis(clr4,gc_cov_4,view_df=view_df_4,n_eigs=3)
cis_eigs_5 = cooltools.eigs_cis(clr1,gc_cov_5,view_df=view_df_5,n_eigs=3)

PC1_timepoint_df = pd.concat([cis_eigs_1[1]['chrom'],cis_eigs_1[1]['start'],cis_eigs_1[1]['end'],cis_eigs_1[1]['E1'],cis_eigs_2[1]['E1'],cis_eigs_3[1]['E1'],cis_eigs_4[1]['E1'],cis_eigs_5[1]['E1']],axis = 1,keys = ['chrom','start','end','control','12hpi','2dpi','4dpi','7dpi'])
PC1_timepoint_df['Chomosome_Region'] = PC1_timepoint_df.apply(lambda row : f"chr{row['chrom']}:{row['start']}-{row['end']}", axis = 1)
PC2_timepoint_df = pd.concat([cis_eigs_1[1]['chrom'],cis_eigs_1[1]['start'],cis_eigs_1[1]['end'],cis_eigs_1[1]['E2'],cis_eigs_2[1]['E2'],cis_eigs_3[1]['E2'],cis_eigs_4[1]['E2'],cis_eigs_5[1]['E2']],axis = 1,keys = ['chrom','start','end','control','12hpi','2dpi','4dpi','7dpi'])
PC2_timepoint_df['Chomosome_Region'] = PC2_timepoint_df.apply(lambda row : f"chr{row['chrom']}:{row['start']}-{row['end']}", axis = 1)
PC3_timepoint_df = pd.concat([cis_eigs_1[1]['chrom'],cis_eigs_1[1]['start'],cis_eigs_1[1]['end'],cis_eigs_1[1]['E3'],cis_eigs_2[1]['E3'],cis_eigs_3[1]['E3'],cis_eigs_4[1]['E3'],cis_eigs_5[1]['E3']],axis = 1,keys = ['chrom','start','end','control','12hpi','2dpi','4dpi','7dpi'])
PC3_timepoint_df['Chomosome_Region'] = PC3_timepoint_df.apply(lambda row : f"chr{row['chrom']}:{row['start']}-{row['end']}", axis = 1)
feather.write_feather(PC1_timepoint_df,f"/home/ksslab/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/Overall_table/Overall_PC1_{resolution}.feather")
feather.write_feather(PC2_timepoint_df,f"/home/ksslab/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/Overall_table/Overall_PC2_{resolution}.feather")
feather.write_feather(PC3_timepoint_df,f"/home/ksslab/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/Overall_table/Overall_PC3_{resolution}.feather")

print(PC1_timepoint_df)
print('\n\n')
print(PC2_timepoint_df)
print('\n\n')
print(PC3_timepoint_df)