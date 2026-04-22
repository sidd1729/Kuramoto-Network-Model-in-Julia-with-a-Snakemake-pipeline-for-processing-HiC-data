#Reading and converting the counts to tpm
#Sys.setenv(http_proxy = '172.16.2.252:3128')
#Sys.setenv(https_proxy = '172.16.2.252:3128')
#install.packages("readxl",repos = "https://cloud.r-project.org")
library(readxl)
library(purrr)
library(furrr)
library(dplyr)
library(arrow)
unique_gene_length_tables <- read.csv("~/Downloads/unique_danRer10_refGene_gene.names.transcript.length.txt",header = TRUE)
data_1 <- read_excel("~/Downloads/GSE180518_Table-1_Kramer_2021.xlsx",sheet = "Pairwise_CtrlTo24hpl")
data_2 <- read_excel("~/Downloads/GSE180518_Table-1_Kramer_2021.xlsx",sheet = "Pairwise_CtrlTo36hpl")
data_3 <- read_excel("~/Downloads/GSE180518_Table-1_Kramer_2021.xlsx",sheet = "Pairwise_CtrlTo72hpl")
data_4 <- read_excel("~/Downloads/GSE180518_Table-1_Kramer_2021.xlsx",sheet = "Pairwise_CtrlTo5dpl")
head(data_1)
head(data_2)
head(data_3)
head(data_4)
data_required <- data_1 %>%
  full_join(data_2, by = "...1") %>%
  full_join(data_3, by = "...1") %>%
  full_join(data_4, by = "...1") %>%
  full_join(unique_gene_length_tables, by = c("...1" = "unique_genes"))

data_subset <- data_required %>%
  select(
    gene_name = ...1,
    gene_lengths,
    starts_with("Ctlr"),
    starts_with("24Hour"),
    starts_with("36Hour"),
    starts_with("72Hour"),
    starts_with("5Day")
  )
write.csv(data_subset,"~/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/RNA_seq_processing/Kramer_2021_subset_raw_counts.csv")

data_normalized <- data_subset %>%
  mutate(across(
    -c(gene_name, gene_lengths), 
    ~ (.x * 1000) / gene_lengths
  ))
write.csv(data_normalized,"~/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/RNA_seq_processing/Kramer_2021_subset_rpk.csv")
data_TPMs <- data_normalized %>% mutate(across(
  -c(gene_name, gene_lengths), 
  ~ (.x * 1000000) / sum(.x,na.rm = TRUE)
))
test_sum <- data_TPMs %>% summarise(across(where(is.numeric),\(x) sum(x,na.rm = TRUE)))
test_var <- data_TPMs %>% summarise(across(where(is.numeric),\(x) sd(x,na.rm = TRUE)))

data_TPMs <- data_TPMs %>%
  mutate(
    ctrl_mean   = rowMeans(across(starts_with("Ctlr")), na.rm = TRUE),
    hpi_24_mean = rowMeans(across(starts_with("24Hour")), na.rm = TRUE),
    hpi_36_mean = rowMeans(across(starts_with("36Hour")), na.rm = TRUE),
    hpi_72_mean = rowMeans(across(starts_with("72Hour")), na.rm = TRUE),
    dpi_5_mean  = rowMeans(across(starts_with("5Day")), na.rm = TRUE)
  )
write.csv(data_TPMs,"~/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/RNA_seq_processing/Kramer_2021_subset_tpm.csv")
processed_df <- data_TPMs[c("gene_name","gene_lengths","ctrl_mean","hpi_24_mean","hpi_36_mean","hpi_72_mean","dpi_5_mean")]
write.csv(processed_df,"~/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/RNA_seq_processing/Kramer_2021_subset_tpm_means.csv")
write_feather(processed_df,"~/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/RNA_seq_processing/Kramer_2021_subset_tpm_means.feather")
gene_name_col <- data_1$...1
ctrl_col_1 <- data_1$Ctlr1
ctrl_col_2 <- data_1$Ctlr2
ctrl_col_3 <- data_1$Ctlr3
ctrl_col_4 <- data_1$Ctlr4
ctrl_col_5 <- data_1$Ctlr5
ctrl_col_6 <- data_1$Ctlr6
hpi_24_col_1 <- data_1$`24Hour1`
hpi_24_col_2 <- data_1$`24Hour2`
hpi_24_col_3 <- data_1$`24Hour3`
hpi_24_col_4 <- data_1$`24Hour4`
hpi_24_col_5 <- data_1$`24Hour5`
hpi_24_col_6 <- data_1$`24Hour6`
hpi_36_col_1 <- data_2$`36Hour1`
hpi_36_col_2 <- data_2$`36Hour2`
hpi_36_col_3 <- data_2$`36Hour3`
hpi_36_col_4 <- data_2$`36Hour4`
hpi_36_col_5 <- data_2$`36Hour5`
hpi_36_col_6 <- data_2$`36Hour6`
hpi_72_col_1 <- data_3$`72Hour1`
hpi_72_col_2 <- data_3$`72Hour2`
hpi_72_col_3 <- data_3$`72Hour3`
hpi_72_col_4 <- data_3$`72Hour4`
hpi_72_col_5 <- data_3$`72Hour5`
hpi_72_col_6 <- data_3$`72Hour6`
dpi_5_col_1 <- data_4$`5Day1`
dpi_5_col_2 <- data_4$`5Day2`
dpi_5_col_3 <- data_4$`5Day3`
dpi_5_col_4 <- data_4$`5Day4`
dpi_5_col_5 <- data_4$`5Day5`
dpi_5_col_6 <- data_4$`5Day6`
data_required <- data.frame(ctr_1 = ctrl_col_1, ctrl_2 = ctrl_col_2, ctrl_3 = ctrl_col_3)

#gene_length_normalizer <- function(i,data_selected = "data_1",column_selected="Ctlr1"){
  #data_selected_1 = get(data_selected)
  #selected_data = data_selected_1[,str_glue("{column_selected}")] 
  #print(selected_data)
  #matched_index <-which(data_selected_1$...1 == unique_gene_length_tables$unique_genes[i])
  #RPM_value = selected_data[matched_index] * 1000 /(unique_gene_length_tables$gene_lengths[i])
  #return(RPM_value)
#}
#colnames(data_1)
#plan(multisession(workers = 30))
#worked <- future_pmap(list(seq(1,dim(data_1)[1],by = 1),rep("data_1",dim(data_1)[1]),rep("Ctlr1",dim(data_1)[1])),gene_length_normalizer)