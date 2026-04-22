#Pipeline for Curve fitting modified with arrow allowing for reading of same data between python and R
library(arrow)
library(purrr) #Used for map function
library(furrr) #Used for parallel map function
library(parallel)
library(dplyr)
library(zoo)
library(stringr) #For string interpolation

#Functions:
K_50_calc_function<-function(K_vals,R_vals,i_val,base_title,base_plot_path,end_plot_string,plot_disp){
  K_vals <- unlist(K_vals)
  R_vals <- unlist(R_vals)
  iteration_no <- 0
  init_up <- 0.06
  if (all(is.na(R_vals))){
    return (list(NaN,NaN))#(list(-1,-1))
  }
  spl <- smooth.spline(na.spline(K_vals),na.spline(R_vals))
  Sampling_seq<-seq(from = 0, to = init_up, length.out = 50000)
  predicted_vals <- predict(spl,Sampling_seq)
  
  desired_index <- which.min(abs(predicted_vals$y-0.5))
  K_50 <- predicted_vals$x[desired_index]
  R_50 <- predicted_vals$y[desired_index]
  

  while ((abs(R_50 - 0.5) > 0.02) && (iteration_no < 1000)) {
    #spl <- smooth.spline(na.spline(K_vals),na.spline(R_vals))
    Sampling_seq<-seq(from = max(0,K_50 - 1/(2^iteration_no)), to = min(init_up,K_50 + 1/(2^iteration_no)), length.out = 5000)
    predicted_vals <- predict(spl,Sampling_seq)
    
    desired_index <- which.min(abs(predicted_vals$y-0.5))
    K_50 <- predicted_vals$x[desired_index]
    R_50 <- predicted_vals$y[desired_index]
    iteration_no <- iteration_no + 1
  }

  while ((abs(R_50 - 0.5) > 0.02) && (iteration_no < 1000)) {
    #spl <- smooth.spline(na.spline(K_vals),na.spline(R_vals))
    Sampling_seq<-seq(from = 0, to = init_up/(2^iteration_no), length.out = 5000)
    predicted_vals <- predict(spl,Sampling_seq)
    
    desired_index <- which.min(abs(predicted_vals$y-0.5))
    K_50 <- predicted_vals$x[desired_index]
    R_50 <- predicted_vals$y[desired_index]
    iteration_no <- iteration_no + 1
  }
  
  if (plot_disp == TRUE){
  invisible(png(str_glue("{base_plot_path}/plot_{i_val}_{end_plot_string}")))
  plot(K_vals,R_vals,pch=19,main = str_glue("{base_title} for index {i_val}"),xlab = 'k',ylab = 'r')
  abline(h = 0.5, col = "blue", lty = 2,lwd = 3)  # Horizontal line
  abline(h = R_50,col = "yellow",lty = 2,lwd = 3)
  abline(v = K_50, col = "red", lty = 3,lwd = 3)   # Vertical line
  text("Predicted",x = K_50 ,y = max(R_vals))
  lines(spl,lwd = 3,col = 'green')
  invisible(dev.off())}
  #cat("\014")
  output_vec <- list(K_50,R_50)
  return (output_vec)}
#mapped <- map(output_vec,print)

Perpetual_1_classifier <- function(K_vals,R_vals,i_val){
  len_improper = length(R_vals[R_vals > 0.5])
  if (len_improper == length(R_vals)){
    return(i_val)
  } else {
    return(-1)
  }
}
#library(tidyr) #To get final output into 2 separate columns
args <- commandArgs(trailingOnly = TRUE)
k <- read_feather(args[1])
r_kura <- read_feather(args[2])
r_uni <- read_feather(args[3])

col_no <- ncol(r_uni)
rows1 <- rep(k,col_no)
rows2 <- r_kura
rows3 <- r_uni

#plan(multisession, workers = as.integer(args[4]),gc=TRUE)
i  = seq_along(rows1)

if (length(args) > 8){
  base_title = args[9]
  base_plot_path = args[10]
  end_plot_string = args[11]
  plot_disp = TRUE
} else {
  base_title = "k-r for ICE-Kura order chr1 100kb 5x5 GM12878"
  base_plot_path = "/home/ksslab/Siddharth/Temp_plot_dump/GM12878_5_5_chr1_kura"
  end_plot_string = ".png"
  plot_disp = FALSE
}
base_title_vec = rep(str_glue({base_title}),length(rows1))
base_plot_vec = rep(str_glue({base_plot_path}),length(rows2))
end_plot_vec = rep(str_glue({end_plot_string}),length(rows2))
plot_disp_vec = rep(plot_disp,length(rows1))

#K_50 calculations
time1<-system.time(Out_1 <-future_pmap(list(rows1,rows2,i,base_title_vec,base_plot_vec,end_plot_vec,plot_disp_vec),K_50_calc_function))
time2<-system.time(Out_2 <-future_pmap(list(rows1,rows3,i,base_title_vec,base_plot_vec,end_plot_vec,plot_disp_vec),K_50_calc_function))
write_feather(tibble(Out_1),args[5])
write_feather(tibble(Out_2),args[6])

#Identifying perpetually 1 graphs
time3<-system.time(Out_perpet_1 <-future_pmap(list(rows1,rows2,i),Perpetual_1_classifier))
time4<-system.time(Out_perpet_2 <-future_pmap(list(rows1,rows3,i),Perpetual_1_classifier))
write_feather(tibbled_out<-tibble(Out_perpet_1),args[7])
write_feather(tibbled_out<-tibble(Out_perpet_2),args[8])


#time1,time2,time3,time4