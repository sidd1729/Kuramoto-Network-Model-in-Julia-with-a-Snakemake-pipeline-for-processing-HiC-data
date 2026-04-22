#Curve manual testing ideally for curve fitting
#library(readr)
#library(psych) #Used for tr function for matrix trace
library(purrr) #Used for map function
library(reticulate)
args <- commandArgs(trailingOnly = TRUE)
Sys.setenv(RETICULATE_PYTHON = "C:/Users/HP/AppData/Local/Programs/Python/Python313/")
#py_require("numpy")    # Declare jax is a requirement
#import_builtins(convert = TRUE, delay_load = FALSE)
#numpy <- import("numpy")
#pickle <- import("pickle")
#np_mat<-numpy$load("C:/Users/HP/Downloads/Siddharth/Data/hg38/stem_cell_H1_HESC/100kb/raw/chr1.mat.txt.npy")
#Remove delete line from matrix_reader for the matrices
#file <-  open('data.pkl', 'rb')
#loaded_data = pickle$load(file)

K_vals_obj<-py_load_object("C:/Users/HP/Downloads/Siddharth/Scripts/Python_scripts/HESC/ICE/10k_coupling_k_vals__500_Chr1_100kb_ICE_GM12878_0.pickle", pickle = "pickle", convert = TRUE)
R_vals_obj<-py_load_object("C:/Users/HP/Downloads/Siddharth/Scripts/Python_scripts/HESC/ICE/10k_coupling_trial_500_Chr1_100kb_ICE_GM12878_0.pickle", pickle = "pickle", convert = TRUE)
R_vals<-unlist(map(R_vals_obj,function(x) x[[1]]))
K_vals<-unlist(map(K_vals_obj,function(x) x[[1]]))

iteration_no <- 0
init_up <- 0.025
spl <- smooth.spline(K_vals,R_vals)
Sampling_seq<-seq(from = 0, to = init_up, length.out = 50000)
predicted_vals <- predict(spl,Sampling_seq)

desired_index <- which.min(abs(predicted_vals$y-0.5))
K_50 <- predicted_vals$x[desired_index]
R_50 <- predicted_vals$y[desired_index]

#while ((abs(R_50 - 0.5) > 0.02) && (iteration_no < 1000)) {
#spl <- smooth.spline(K_vals,R_vals)
#Sampling_seq<-seq(from = 0, to = init_up/(2^iteration_no), length.out = 50000)
#predicted_vals <- predict(spl,Sampling_seq)

#desired_index <- which.min(abs(predicted_vals$y-0.5))
#K_50 <- predicted_vals$x[desired_index]
#R_50 <- predicted_vals$y[desired_index]
#iteration_no <- iteration_no + 1
#}

#svg("C:/Users/HP/Downloads/plot.svg")
invisible(png("C:/Users/HP/Downloads/plot_0_GM12878_ICE.png"))
plot(K_vals,R_vals,pch=19,main = "k-r for Raw-Kura order chr1 100kb index 0 of 10x10 GM12878 ICE", xlab='k', ylab='r')
abline(h = 0.5, col = "blue", lty = 2,lwd = 3)  # Horizontal line
abline(h = R_50,col = "yellow",lty = 2,lwd = 3)
abline(v = K_50, col = "red", lty = 3,lwd = 3)   # Vertical line
text("Predicted",x = K_50 ,y = max(R_vals))
lines(spl,lwd = 3,col = 'green')
invisible(dev.off())
#cat("\014")
output_vec <- c(K_50,R_50)
mapped <- map(output_vec,print)
