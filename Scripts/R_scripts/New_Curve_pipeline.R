#Pipeline for Curve fitting modified with arrow allowing for reading of same data between python and R
library(arrow)
library(purrr) #Used for map function
library(furrr) #Used for parallel map function
library(parallel)
library(dplyr)
library(zoo)
library(stringr) #For string interpolation

#Functions:
K_50_calc_function<-function(K_vals,R_vals,i_val,base_title,base_plot_path,end_plot_string){
  K_vals <- unlist(K_vals)
  R_vals <- unlist(R_vals)
  iteration_no <- 0
  init_up <- 0.025
  if (all(is.na(R_vals))){
    return (list(-1,-1))#(list(NA,NA))
  }
  spl <- smooth.spline(na.spline(K_vals),na.spline(R_vals))
  Sampling_seq<-seq(from = 0, to = init_up, length.out = 50000)
  predicted_vals <- predict(spl,Sampling_seq)
  
  desired_index <- which.min(abs(predicted_vals$y-0.5))
  K_50 <- predicted_vals$x[desired_index]
  R_50 <- predicted_vals$y[desired_index]
  
  while ((abs(R_50 - 0.5) > 0.02) && (iteration_no < 1000)) {
    spl <- smooth.spline(na.spline(K_vals),na.spline(R_vals))
    Sampling_seq<-seq(from = 0, to = init_up/(2^iteration_no), length.out = 50000)
    predicted_vals <- predict(spl,Sampling_seq)
    
    desired_index <- which.min(abs(predicted_vals$y-0.5))
    K_50 <- predicted_vals$x[desired_index]
    R_50 <- predicted_vals$y[desired_index]
    iteration_no <- iteration_no + 1
  }
  
  #svg("C:/Users/HP/Downloads/plot.svg")
  #invisible(png(args[4]))
  #invisible(png(str_glue("{base_plot_path}/plot_{i_val}_GM12878_ICE_5_5_Kura.png"))) #Main line
  #invisible(png(str_glue("{base_plot_path}/plot_{i_val}_{end_plot_string}))) #Main line 2
  plot(K_vals,R_vals,pch=19,main = str_glue("{base_title} for index {i_val}"),xlab = 'k',ylab = 'r')#args[5], xlab='k', ylab='r')
  abline(h = 0.5, col = "blue", lty = 2,lwd = 3)  # Horizontal line
  abline(h = R_50,col = "yellow",lty = 2,lwd = 3)
  abline(v = K_50, col = "red", lty = 3,lwd = 3)   # Vertical line
  text("Predicted",x = K_50 ,y = max(R_vals))
  lines(spl,lwd = 3,col = 'green')
  invisible(dev.off())
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
feathered <- read_feather(args[1])#("/home/ksslab/Siddharth/Program_Directory/Data/Simulation_Results/Linux/GM12878/143_158_Reproducibility/Chr1_Kura_k_1.feather")
feathered_r <- read_feather(args[2])#("/home/ksslab/Siddharth/Program_Directory/Data/Simulation_Results/Linux/GM12878/143_158_Reproducibility/Chr1_Kura_r_1.feather")
feathered_comb <- cbind(feathered,feathered_r)
rm(feathered)
rm(feathered_r)
col_no <- ncol(feathered_comb)
colnames(feathered_comb) <- 1:col_no
#feathered_comb_transp <- transpose(feathered_comb)
feathered_comb<- feathered_comb %>% rowwise() %>% mutate(across(everything(),function(x) x[[1]]))
rows1 <- split(feathered_comb[,1: (col_no %/%2) ], seq(nrow(feathered_comb[,1: (col_no %/%2) ])))
rows2 <- split(feathered_comb[,((col_no%/%2) + 1):col_no], seq(nrow(feathered_comb[,((col_no%/%2) + 1):col_no])))
#mutate(across(everything(),function(x) x[[1]]))
#map2(feathered_comb[,1: (col_no %/%2) ] %>% rowwise() ,feathered_comb[,((col_no%/%2) + 1):col_no] %>% rowwise(),sum)
rm(feathered_comb)
rows2[1:2]
#plan(multisession, workers = 30)
#time1<-system.time(Out <-future_map2(rows1[1:10],rows2[1:10],K_50_calc_function))
plan(multisession, workers = 30,gc=TRUE)
i  = seq_along(rows1)
length(rows1)
base_title = args[3]#"k-r for ICE-Kura order chr1 100kb 5x5 GM12878"
base_title_vec = rep(str_glue({base_title}),length(rows1))
base_plot_path = args[4]#"/home/ksslab/Siddharth/Temp_plot_dump/GM12878_5_5_chr1_kura"
base_plot_vec = rep(str_glue({base_plot_path}),length(rows2))
time2<-system.time(Out <-future_pmap(list(rows1,rows2,i,base_title_vec,base_plot_vec),K_50_calc_function))
write_feather(tibble(Out),args[5])#"/home/ksslab/Siddharth/Program_Directory/Data/Simulation_Results/Linux/GM12878/143_158_Reproducibility/Chr1_K_50_Kura_1.feather")


time1<-system.time(Out_perpet_1 <-future_pmap(list(rows1,rows2,i),Perpetual_1_classifier))
write_feather(tibbled_out<-tibble(Out_perpet_1),args[6])#"/home/ksslab/Siddharth/Program_Directory/Data/Simulation_Results/Linux/GM12878/143_158_Reproducibility/Chr1_perpet_1_Kura_1.feather")

end_plot_string = args[7]
#Previous was tibble(Out)
#work_pls <-tibble(data=Out)
#out_of_a_bind <- do.call(rbind, Out)
#out_tibbled <- tibble(out_of_a_bind)
#data.table::rbindlist(l = lapply(Out, data.table::rbindlist, fill = TRUE), fill = TRUE)
#unnest(tibble(Out))
#length(rows1[1:2])
#feathered$`0`
#print(feathered_comb %>% rowwise() %>% select(1:500))
#cl <- makeCluster(2)
#time_1<-system.time(
#trialed <- apply(feathered,MARGIN = 1,function(x) x[[1]]))

#time_2 <- system.time(trialed_parallel <- parApply(cl,feathered,MARGIN = 1,function(x) x[[1]]))

#Curve_pipeline_fn(){
  
#}
#library(parallelly)
#if (supportsMulticore() == TRUE) {
  #print("Yay")
#}  else {
  #print("Nay")
#}


