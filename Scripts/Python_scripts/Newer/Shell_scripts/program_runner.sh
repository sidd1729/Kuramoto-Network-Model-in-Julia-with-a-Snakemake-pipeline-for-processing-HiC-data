#!/bin/bash

# Control S63
#echo "Control S63"
#python Repeated_recalc_K_50.py /home/ksslab/Siddharth/Program_Directory/Data/Intermediate_Processing_Data/control_S63.feather /home/ksslab/Siddharth/Program_Directory/Data/Intermediate_Processing_Data/S63_K_50/Kura_r /home/ksslab/Siddharth/Program_Directory/Data/Intermediate_Processing_Data/S63_K_50/Kura_k /home/ksslab/Siddharth/Program_Directory/Data/Intermediate_Processing_Data/S63_K_50/Uni_r /home/ksslab/Siddharth/Program_Directory/Data/Intermediate_Processing_Data/S63_K_50/Uni_k
#echo "Control S63 finished."

#sleep 3h

# Control S68
#echo "Starting Control S68"
#python Repeated_recalc_K_50.py /home/ksslab/Siddharth/Program_Directory/Data/Intermediate_Processing_Data/control_S68.feather /home/ksslab/Siddharth/Program_Directory/Data/Intermediate_Processing_Data/S_68_K_50/Kura_r /home/ksslab/Siddharth/Program_Directory/Data/Intermediate_Processing_Data/S_68_K_50/Kura_k /home/ksslab/Siddharth/Program_Directory/Data/Intermediate_Processing_Data/S_68_K_50/Uni_r /home/ksslab/Siddharth/Program_Directory/Data/Intermediate_Processing_Data/S_68_K_50/Uni_k
#echo "Control S68 finished."

#sleep 3h

# 12_hpi_ S165
#echo "Starting 12 hpi S165"
#python Repeated_recalc_K_50.py /home/ksslab/Siddharth/Program_Directory/Data/Intermediate_Processing_Data/12_hpi_S165.feather /home/ksslab/Siddharth/Program_Directory/Data/Intermediate_Processing_Data/S165_K_50/Kura_r /home/ksslab/Siddharth/Program_Directory/Data/Intermediate_Processing_Data/S165_K_50/Kura_k /home/ksslab/Siddharth/Program_Directory/Data/Intermediate_Processing_Data/S165_K_50/Uni_r /home/ksslab/Siddharth/Program_Directory/Data/Intermediate_Processing_Data/S165_K_50/Uni_k
#echo "12 hpi S165 completed"

#sleep 3h

# 12 hpi S207
#echo "Starting 12 hpi S207"
#python Repeated_recalc_K_50.py /home/ksslab/Siddharth/Program_Directory/Data/Intermediate_Processing_Data/12_hpi_S207.feather /home/ksslab/Siddharth/Program_Directory/Data/Intermediate_Processing_Data/S207_K_50/Kura_r /home/ksslab/Siddharth/Program_Directory/Data/Intermediate_Processing_Data/S207_K_50/Kura_k /home/ksslab/Siddharth/Program_Directory/Data/Intermediate_Processing_Data/S207_K_50/Uni_r /home/ksslab/Siddharth/Program_Directory/Data/Intermediate_Processing_Data/S207_K_50/Uni_k
#echo "12 hpi S207 completed"

#sleep 3h

#echo "All programs completed."

echo "Starting control 63"
Rscript ~/Siddharth/Program_Directory/Scripts/R_scripts/New_Curve_pipeline.R /home/ksslab/Siddharth/Program_Directory/Data/Intermediate_Processing_Data/S63_K_50/Kura_k /home/ksslab/Siddharth/Program_Directory/Data/Intermediate_Processing_Data/S63_K_50/Kura_r "k-r for ICE-Kura order chr1 40kb 5x5 15Mb control 63" /home/ksslab/Siddharth/Temp_plot_dump/GM12878_5_5_chr1_kura /home/ksslab/Siddharth/Program_Directory/Data/Intermediate_Processing_Data/S63_K_50/S63_K50 /home/ksslab/Siddharth/Program_Directory/Data/Intermediate_Processing_Data/S63_K_50/Perpet_1
echo "Control S63 completed"

echo "Starting control 68"
Rscript ~/Siddharth/Program_Directory/Scripts/R_scripts/New_Curve_pipeline.R /home/ksslab/Siddharth/Program_Directory/Data/Intermediate_Processing_Data/S_68_K_50/Kura_k /home/ksslab/Siddharth/Program_Directory/Data/Intermediate_Processing_Data/S_68_K_50/Kura_r "k-r for ICE-Kura order chr1 40kb 5x5 15Mb control 68" /home/ksslab/Siddharth/Temp_plot_dump/GM12878_5_5_chr1_kura /home/ksslab/Siddharth/Program_Directory/Data/Intermediate_Processing_Data/S_68_K_50/S68_K50 /home/ksslab/Siddharth/Program_Directory/Data/Intermediate_Processing_Data/S_68_K_50/Perpet_1
echo "Control S68 completed"

echo "Starting 12hpi S165"
Rscript ~/Siddharth/Program_Directory/Scripts/R_scripts/New_Curve_pipeline.R /home/ksslab/Siddharth/Program_Directory/Data/Intermediate_Processing_Data/S165_K_50/Kura_k /home/ksslab/Siddharth/Program_Directory/Data/Intermediate_Processing_Data/S165_K_50/Kura_r "k-r for ICE-Kura order chr1 40kb 5x5 15Mb 12hpi S165" /home/ksslab/Siddharth/Temp_plot_dump/GM12878_5_5_chr1_kura /home/ksslab/Siddharth/Program_Directory/Data/Intermediate_Processing_Data/S165_K_50/S165_K50 /home/ksslab/Siddharth/Program_Directory/Data/Intermediate_Processing_Data/S165_K_50/Perpet_1
echo "12 hpi S165 completed"

echo "Starting 12hpi S207"
Rscript ~/Siddharth/Program_Directory/Scripts/R_scripts/New_Curve_pipeline.R /home/ksslab/Siddharth/Program_Directory/Data/Intermediate_Processing_Data/S207_K_50/Kura_k /home/ksslab/Siddharth/Program_Directory/Data/Intermediate_Processing_Data/S207_K_50/Kura_r "k-r for ICE-Kura order chr1 40kb 5x5 15Mb 12hpi S207" /home/ksslab/Siddharth/Temp_plot_dump/GM12878_5_5_chr1_kura /home/ksslab/Siddharth/Program_Directory/Data/Intermediate_Processing_Data/S207_K_50/S207_K50 /home/ksslab/Siddharth/Program_Directory/Data/Intermediate_Processing_Data/S207_K_50/Perpet_1
echo "12 hpi S207 completed"