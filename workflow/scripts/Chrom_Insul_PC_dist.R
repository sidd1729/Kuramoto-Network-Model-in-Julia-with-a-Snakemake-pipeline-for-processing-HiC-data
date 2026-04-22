library(arrow)
library(stringr)
resolution = 50000
count = 0
chromosomes <- seq(1,25,by=1)
corres_vector = c()
colors_vec <- colors()[seq(50,500,by=20)]
colors_vec <- c("blue","green","orange","red","purple","blue","yellow")
for (chromosome in chromosomes){
count = 0
corres_vector = c()
colors_vec <- colors()[seq(50,500,by=20)]
colors_vec <- c("blue","green","orange","red","purple","blue","yellow")
svg(filename = str_glue("/home/ksslab/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/temporaries/trials/R_trials/{resolution}/{chromosome}/PC1_{resolution}_{chromosome}.svg"))#,res = 300)
#png(filename = str_glue("/home/ksslab/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/temporaries/trials/R_trials/{resolution}/{chromosome}/PC3_{resolution}_{chromosome}.png"),res = 200)
#svg(str_glue("/home/ksslab/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/temporaries/trials/R_trials/{resolution}/Rao_{resolution}.svg"))s
#file = str_glue("/home/ksslab/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/temporaries/trials/4DNFITRVKRPA_{resolution}_PC1.feather")
for (file in list.files(str_glue("/home/ksslab/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/PCs/{resolution}"))){
    #print(file)
    #print('\n\n')
    #print(str_split(file,'_'))
    #print(str_split(file,'_')[[1]])
    #print(str_split(file,'_')[[1]][length(str_split(file,'_')[[1]])-1])
    #print(substr(file, str_length(file)-10, str_length(file)-8))
    if ((substr(file, str_length(file)-10, str_length(file)-8) ==  "PC1") & (str_split(file,'_')[[1]][length(str_split(file,'_')[[1]])-1] == str_glue('{chromosome}')) & (str_split(file,'_')[[1]][length(str_split(file,'_')[[1]])-2] == str_glue('{resolution}'))){
    #if (file[-11:-8:+1] == "PC1") & 
        #print(file)
        count <- count + 1
        trial <- read_feather(str_glue("/home/ksslab/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/PCs/{resolution}/{file}"))
        #trial <- read_feather(file)
        trial_nums <- trial$E1 #Change to reflect number of PCs
        trial_nanless <- trial_nums[!is.na(trial_nums)]
        density_plot <- density(trial_nanless)
        if (count == 1){
        #par(mar = c(2,2,2,2))
        plot(density_plot$x,density_plot$y,type = "l",col = colors_vec[count],xlab ="PC1",ylab = "Density of PC1",ylim = c(0,2.5),cex.axis = 0.5)}
        else{
            lines(density_plot$x,density_plot$y,col = colors_vec[count])
        }
         corres_vector<-append(corres_vector,str_glue("{str_split(file,'_')[[1]][1]}_{str_split(file,'_')[[1]][2]}:{colors_vec[count]}"))
         #corres_vector<-append(corres_vector,str_glue("{str_split(file,'_')[[1]][1]}_{str_split(file,'_')[[1]][2]}:{colors_vec[count]}"))
        }
}

title(str_glue("PC1 distribution for various time points for chromosome {chromosome} at {resolution} bin size"))
legend("topright",legend = corres_vector,col = colors_vec[seq(1,count,by=1)])
invisible(dev.off())
}