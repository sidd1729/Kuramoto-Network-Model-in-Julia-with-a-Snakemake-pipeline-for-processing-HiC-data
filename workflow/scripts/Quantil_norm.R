library(arrow)
library(stringr)
library(preprocessCore)
library(dplyr)

#chromosome_no = 1
chromosomes <- seq(1,25,by=1)
for (iced in c(TRUE,FALSE)){
    for (chromosome_no in chromosomes){
if (iced == TRUE){
    path_var = "ice_feather"
    path_end = "_ICEd"
    name_var = "iced"
}
else{
    path_var = "raw_feather"
    path_end = ""
    name_var = "raw"
}
trial_1 <- as.matrix(read_feather(str_glue("/home/ksslab/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/temporaries/{path_var}/50000/control_S68_50000_{chromosome_no}{path_end}.feather")))
trial_shape_1 = dim(trial_1)
trial_vec_1 <- as.vector(trial_1)

trial_2 <- as.matrix(read_feather(str_glue("/home/ksslab/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/temporaries/{path_var}/50000/12_hpi_S207_50000_{chromosome_no}{path_end}.feather")))
trial_shape_2 = dim(trial_2)
trial_vec_2 <- as.vector(trial_2)

trial_3 <- as.matrix(read_feather(str_glue("/home/ksslab/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/temporaries/{path_var}/50000/2_dpi_S208_50000_{chromosome_no}{path_end}.feather")))
trial_shape_3 = dim(trial_3)
trial_vec_3 <- as.vector(trial_3)

trial_4 <- as.matrix(read_feather(str_glue("/home/ksslab/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/temporaries/{path_var}/50000/4_dpi_S69_50000_{chromosome_no}{path_end}.feather")))
trial_shape_4 = dim(trial_4)
trial_vec_4 <- as.vector(trial_4)

trial_5 <- as.matrix(read_feather(str_glue("/home/ksslab/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/temporaries/{path_var}/50000/7_dpi_S1_50000_{chromosome_no}{path_end}.feather")))
trial_shape_5 = dim(trial_5)
trial_vec_5 <- as.vector(trial_5)

mat_mats<-do.call(cbind,list(trial_vec_1,trial_vec_2,trial_vec_3,trial_vec_4,trial_vec_5))
norm_quant <- normalize.quantiles(mat_mats)
print(trial_shape_2)
print('\n')
print(length(trial_vec_2))
print(dim(mat_mats))

norm_quant_1 <- norm_quant[,1]
norm_quant_2 <- norm_quant[,2]
norm_quant_3 <- norm_quant[,3]
norm_quant_4 <- norm_quant[,4]
norm_quant_5 <- norm_quant[,5]
print(list(length(norm_quant_1),length(norm_quant_2),length(norm_quant_3),length(norm_quant_4),length(norm_quant_5)))

norm_1 <- matrix(norm_quant_1,nrow = trial_shape_1[1],ncol = trial_shape_1[2])
norm_2 <- matrix(norm_quant_2,nrow = trial_shape_2[1],ncol = trial_shape_2[2])
norm_3 <- matrix(norm_quant_3,nrow = trial_shape_3[1],ncol = trial_shape_3[2])
norm_4 <- matrix(norm_quant_4,nrow = trial_shape_4[1],ncol = trial_shape_4[2])
norm_5 <- matrix(norm_quant_5,nrow = trial_shape_5[1],ncol = trial_shape_5[2])

write_feather(as.data.frame(trial_1),str_glue("/home/ksslab/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/temporaries/trials/Quantile/{chromosome_no}/{name_var}/unquant_{chromosome_no}_1_{name_var}.feather"))
write_feather(as.data.frame(norm_1),str_glue("/home/ksslab/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/temporaries/trials/Quantile/{chromosome_no}/{name_var}/quant_{chromosome_no}_1_{name_var}.feather"))
write_feather(as.data.frame(trial_2),str_glue("/home/ksslab/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/temporaries/trials/Quantile/{chromosome_no}/{name_var}/unquant_{chromosome_no}_2_{name_var}.feather"))
write_feather(as.data.frame(norm_2),str_glue("/home/ksslab/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/temporaries/trials/Quantile/{chromosome_no}/{name_var}/quant_{chromosome_no}_2_{name_var}.feather"))
write_feather(as.data.frame(trial_3),str_glue("/home/ksslab/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/temporaries/trials/Quantile/{chromosome_no}/{name_var}/unquant_{chromosome_no}_3_{name_var}.feather"))
write_feather(as.data.frame(norm_3),str_glue("/home/ksslab/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/temporaries/trials/Quantile/{chromosome_no}/{name_var}/quant_{chromosome_no}_3_{name_var}.feather"))
write_feather(as.data.frame(trial_4),str_glue("/home/ksslab/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/temporaries/trials/Quantile/{chromosome_no}/{name_var}/unquant_{chromosome_no}_4_{name_var}.feather"))
write_feather(as.data.frame(norm_4),str_glue("/home/ksslab/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/temporaries/trials/Quantile/{chromosome_no}/{name_var}/quant_{chromosome_no}_4_{name_var}.feather"))
write_feather(as.data.frame(trial_5),str_glue("/home/ksslab/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/temporaries/trials/Quantile/{chromosome_no}/{name_var}/unquant_{chromosome_no}_5_{name_var}.feather"))
write_feather(as.data.frame(norm_5),str_glue("/home/ksslab/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/temporaries/trials/Quantile/{chromosome_no}/{name_var}/quant_{chromosome_no}_5_{name_var}.feather"))
    }
}
