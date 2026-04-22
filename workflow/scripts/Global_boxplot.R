library(ggplot2)
library(arrow)
library(dplyr)
library(rlang)
sample_wise_df <- read_feather("results/k_50_table/k_50_sample_wise_df.feather")
df_to_hmap <- filter(sample_wise_df,!!sym(names(sample_wise_df)[1]) != -1 & !!sym(names(sample_wise_df)[2]) != -1 & !!sym(names(sample_wise_df)[3]) != -1 & !!sym(names(sample_wise_df)[4]) != -1 & !!sym(names(sample_wise_df)[5]) != -1)
png("results/k_50_table/k_50_wo_perpet_1_boxplot.png")
ggplot(stack(df_to_hmap), aes(x = factor(ind, levels = names(df_to_hmap)), y = values, fill = ind)) + geom_boxplot() +
  labs(title="Distribution of K_50 across timepoints ",x="Sample", y = "K_50 values") + theme(plot.title = element_text(size = 20))
invisible(dev.off())

sample_wise_df <- read_feather("results/k_50_table/k_50_sample_wise_df.feather")
sample_wise_df[sample_wise_df == -1] <- NA
df_to_hmap <- na.omit(sample_wise_df)
png("results/k_50_table/k_50_wo_perpet_1_na_boxplot.png")
ggplot(stack(df_to_hmap), aes(x = factor(ind, levels = names(df_to_hmap)), y = values, fill = ind)) + geom_boxplot() +
  labs(title="Distribution of K_50 across timepoints ",x="Sample", y = "K_50 values") + theme(plot.title = element_text(size = 20))
invisible(dev.off())

insulation_df <- read_feather("results/Overall_table/Overall_insulation.feather")
df_int <- na.omit(insulation_df) 
df_to_hmap <- df_int[names(df_int)[c(4,5,6,7,8)]]
png("results/Overall_table_plot/insulation_boxplot.png")
ggplot(stack(df_to_hmap), aes(x = factor(ind, levels = names(df_to_hmap)), y = values, fill = ind)) + geom_boxplot() +
  labs(title="Distribution of Insulation across timepoints ",x="Sample", y = "Insulation values") + theme(plot.title = element_text(size = 20))
invisible(dev.off())

DLR_df <- read_feather("results/Overall_table/Overall_DLR.feather")
df_int <- na.omit(DLR_df)
df_to_hmap <- df_int[names(df_int)[c(4,5,6,7,8)]]
png("results/Overall_table_plot/DLR_boxplot.png")
ggplot(stack(df_to_hmap), aes(x = factor(ind, levels = names(df_to_hmap)), y = values, fill = ind)) + geom_boxplot() +
  labs(title="Distribution of DLR across timepoints ",x="Sample", y = "DLR values") + theme(plot.title = element_text(size = 20))
invisible(dev.off())

PC1_df <- read_feather("results/Overall_table/Overall_PC1.feather")
df_int <- na.omit(PC1_df)
df_to_hmap <- df_int[names(df_int)[c(4,5,6,7,8)]]
png("results/Overall_table_plot/PC1_boxplot.png")
ggplot(stack(df_to_hmap), aes(x = factor(ind, levels = names(df_to_hmap)), y = values, fill = ind)) + geom_boxplot() +
  labs(title="Distribution of PC1 across timepoints ",x="Sample", y = "PC1 values") + theme(plot.title = element_text(size = 20))
invisible(dev.off())

PC2_df <- read_feather("results/Overall_table/Overall_PC2.feather")
df_int <- na.omit(PC2_df)
df_to_hmap <- df_int[names(df_int)[c(4,5,6,7,8)]]
png("results/Overall_table_plot/PC2_boxplot.png")
ggplot(stack(df_to_hmap), aes(x = factor(ind, levels = names(df_to_hmap)), y = values, fill = ind)) + geom_boxplot() +
  labs(title="Distribution of PC2 across timepoints ",x="Sample", y = "PC2 values") + theme(plot.title = element_text(size = 20))
invisible(dev.off())

PC3_df <- read_feather("results/Overall_table/Overall_PC3.feather")
df_int <- na.omit(PC3_df)
df_to_hmap <- df_int[names(df_int)[c(4,5,6,7,8)]]
png("results/Overall_table_plot/PC3_boxplot.png")
ggplot(stack(df_to_hmap), aes(x = factor(ind, levels = names(df_to_hmap)), y = values, fill = ind)) + geom_boxplot() +
  labs(title="Distribution of PC3 across timepoints ",x="Sample", y = "PC3 values") + theme(plot.title = element_text(size = 20))
invisible(dev.off())