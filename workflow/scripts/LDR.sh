#cd /home/ksslab/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/temporaries/hics
##echo ${pwd}
#for file in *.hic;do
##"/home/ksslab/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/temporaries/hics/2_dpi_S208.hic"
##file=2_dpi_S208.hic;
##chrno=25;
#filebase=$(basename $file .hic); 
##echo "$file";
    #for chrno in {1..25};do
        #h1d one DLR "/home/ksslab/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/temporaries/hics/$file" 50000 chr${chrno} --maxchr 25 --o /home/ksslab/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/LDR/${filebase}_${chrno}_ldr --datatype rawhic --gt /home/ksslab/Siddharth/Program_Directory/Data/Actual_Data/Chromsizes/danRer11.chrom.sizes;
    #done
#done

##
#for file in *.hic;do
#"/home/ksslab/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/temporaries/hics/2_dpi_S208.hic"
#file=2_dpi_S208.hic;
#chrno=25;
#filebase=$(basename $file .hic); 
#echo "$file";
    #for chrno in {1..25};do
        #h1d one DLR "/home/ksslab/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/temporaries/hics/$file" 50000 chr${chrno} --maxchr 25 --o /home/ksslab/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/LDR/${filebase}_${chrno}_ldr --datatype rawhic --gt /home/ksslab/Siddharth/Program_Directory/Data/Actual_Data/Chromsizes/danRer11.chrom.sizes;
    #done
#done

cd /home/ksslab/Downloads
file=4DNFITRVKRPA.hic
#filebase=$(basename $file .hic); 
#for chrno in {1..22};do$
#h1d one DLR "/home/ksslab/Downloads/chr1.txt.gz" 50000 chr1 -o /home/ksslab/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/temporaries/trials/GM12878/1/results/ldr/chr1_ldr --datatype matrix #--gt /home/ksslab/Siddharth/Program_Directory/Data/hg38/hg38.chrom.sizes;
h1d one DLR "/home/ksslab/Downloads/4DNFI9YAVTI1.hic" 50000 chr1 -o /home/ksslab/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/temporaries/trials/GM12878/1/results/ldr/chr1_ldr --datatype rawhic --gt /home/ksslab/Siddharth/Program_Directory/Data/hg38/hg38.chrom.sizes;
#done