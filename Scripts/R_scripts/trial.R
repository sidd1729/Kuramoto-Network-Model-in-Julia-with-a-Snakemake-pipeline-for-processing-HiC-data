#https://outdated.legendu.net/blog/install-r-packages-on-windows-behind-proxy/
#https://bioconductor.org/packages/release/bioc/manuals/megadepth/man/megadepth.pdf
Sys.setenv(http_proxy = '172.16.2.252:3128')
Sys.setenv(https_proxy = '172.16.2.252:3128')
library("megadepth")
example_bw <- system.file("tests", "test.bam.all.bw",
                          package = "megadepth", mustWork = TRUE
)
annotation_file <- system.file("tests", "testbw2.bed",
                               package = "megadepth", mustWork = TRUE
)
Actual_bw <- system.file("C:/User/HP/Downloads","ENCFF691LFQ.bigWig","packa")
bw_cov <- get_coverage("C:/User/HP/Downloads/ENCFF691LFQ.bigWig")#, op = "mean")#, annotation = "C:/Users/HP/Downloads/ENCFF153VOQ.bed")

bed_trial <- read.table("C:/Users/HP/Downloads/ENCFF153VOQ.bed",header = FALSE)
bed_3 <-bed_trial[,1:3]
write.table(bed_3,"C:/Users/HP/Downloads/ENCFF153VOQ_trimmed.bed",sep = "\t",quote = FALSE,row.names = FALSE,col.names = FALSE)
bw_cov <- get_coverage("C:/User/HP/Downloads/ENCFF691LFQ.bigWig", op = "mean", annotation = "C:/Users/HP/Downloads/ENCFF153VOQ_trimmed.bed")
bw_summ <- get_summary("C:/User/HP/Downloads/ENCFF691LFQ.bigWig", op = "mean", annotation = "C:/Users/HP/Downloads/ENCFF153VOQ_trimmed.bed")

#library(purrr)
#library(dynutils)
#library(arrow)
#li <- list(
  #list(1, 1, 1),
  #list(2, 1, 1)
#)
library(arrow)
library(purrr) #Used for map function
library(furrr) #Used for parallel map function
library(parallel)
library(dplyr)
library(zoo)
library(stringr) #For string interpolation
r_uni <- read_feather("/home/ksslab/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/temporaries/r_uni.feather")
r_kura <- read_feather("/home/ksslab/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/temporaries/r_kura.feather")
k <- read_feather("/home/ksslab/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/temporaries/k.feather")

feathered_comb <- cbind(k,r_kura)
rm(feathered)
rm(feathered_r)
col_no <- ncol(r_uni)
colnames(feathered_comb) <- 1:col_no
#feathered_comb_transp <- transpose(feathered_comb)

rows1 <- rep(k,col_no)
rows2 <- r_kura
#mutate(across(everything(),function(x) x[[1]]))
#map2(feathered_comb[,1: (col_no %/%2) ] %>% rowwise() ,feathered_comb[,((col_no%/%2) + 1):col_no] %>% rowwise(),sum)
rm(feathered_comb)
rows2[1:2]
#plan(multisession, workers = 30)
#time1<-system.time(Out <-future_map2(rows1[1:10],rows2[1:10],K_50_calc_function))
plan(multisession, workers = 30,gc=TRUE)
i  = seq_along(rows1)
length(rows1)
base_title = "k-r for ICE-Kura order chr1 100kb 5x5 GM12878"
base_title_vec = rep(str_glue({base_title}),length(rows1))
base_plot_path = "/home/ksslab/Siddharth/Temp_plot_dump/GM12878_5_5_chr1_kura"
base_plot_vec = rep(str_glue({base_plot_path}),length(rows2))

K_50_calc_function<-function(K_vals,R_vals,i_val,base_title,base_plot_path,end_plot_string){
  K_vals <- unlist(K_vals)
  #print(K_vals)
  R_vals <- unlist(R_vals)
  #print(R_vals)
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
time2<-system.time(Out <-future_pmap(list(rows1,rows2,i,base_title_vec,base_plot_vec),K_50_calc_function))
#tib <- list_as_tibble(li)

unlist(rows1$data)
