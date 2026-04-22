

#conda init zsh
#conda activate cooltools

cd /home/ksslab/Siddharth/Program_Directory/Data/Actual_Data/Raw_pairs

for file in *.bsorted.pairs; 
do filebase=$(basename $file .bsorted.pairs); 
echo ${filebase}; 
java -jar /home/ksslab/Programs/Juicer/juicer_tools_1.22.01.jar pre $file /home/ksslab/Siddharth/Program_Directory/Data/Actual_Data/Raw_hic_corrected/${filebase}.hic /home/ksslab/Siddharth/Program_Directory/Data/Actual_Data/Chromsizes/DanRer11_orig_HiCup_chrom.size.txt
hic2cool convert /home/ksslab/Siddharth/Program_Directory/Data/Actual_Data/Raw_hic_corrected/${filebase}.hic /home/ksslab/Siddharth/Program_Directory/Data/Actual_Data/Raw_cool_corrected/${filebase}.mcool
done





