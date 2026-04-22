#Yet another density script
library(arrow)
library(dplyr)
library(stringr)
chromosome = 19
chr_selected = str_glue("chr{chromosome}")
resolution = "50kb"
tester <- read_feather("~/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/Cooltools_compartment_cluster/dcHiC_Regen_filtered_non_clustered_compartments_yan.feather")
print(head(tester))
tester_chr_x <- tester %>% filter(chr == str_glue("{chr_selected}"))
print(head(tester_chr_x))
control_den <- density(tester_chr_x$control_100kb)
hpi_12_den <- density(tester_chr_x$hpi_12_100kb)
dpi_2_den <- density(tester_chr_x$dpi_2_100kb)
dpi_4_den <- density(tester_chr_x$dpi_4_100kb)
dpi_7_den <- density(tester_chr_x$dpi_7_100kb)

#Making density plots
svg(str_glue("~/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/Cooltools_compartment_cluster/dcHiC_{chr_selected}_density.svg"))#,width = 20, height = 20)
plot(control_den$x,control_den$y,type = "l",col = "purple",xlab = str_glue("{chr_selected} dcHiC PC"),ylab = "Density of dcHiC PC",ylim = c(0,0.06),cex.axis = 0.5)
lines(hpi_12_den$x,hpi_12_den$y,col = "blue")
lines(dpi_2_den$x,dpi_2_den$y,col = "green")
lines(dpi_4_den$x,dpi_4_den$y,col = "orange")
lines(dpi_7_den$x,dpi_7_den$y,col = "red")
title(main = str_glue("Distribution of dcHiC for various time points for chromosome {chromosome} at {resolution} bin size"),cex.main = 0.8)
legend("topright",legend = c("control","12hpi","2dpi","4dpi","7dpi"),col = c("purple","blue","green","orange","red"),lty = 19)
dev.off()