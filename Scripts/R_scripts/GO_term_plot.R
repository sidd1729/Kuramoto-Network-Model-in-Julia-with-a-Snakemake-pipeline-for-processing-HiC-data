library(ggplot2)
library(stringr)
library(dplyr)
library(purrr)
nature <- "raw"
cluster_max <- 8 #set to 8 when using the mfuzz cluster 32 when sign clusterm
type <- "dc_thres" #"cool" "two_tp_tad"#
thres_decided = 0.2592
two_timepoint_options <- c("control_hpi_12","hpi_12_dpi_2","dpi_2_dpi_4","dpi_4_dpi_7")
two_timepoint_selected <- two_timepoint_options[2]

cluster_plotter <- function(cluster){
if(type == "cool"){
go_readout <- read.csv(str_glue("~/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/Cooltools_compartment_cluster/clustered_cooltools_PCs_sample_wise_{nature}_k_means_cluster_{cluster}_GO_BP.csv"))
} else if (type == "dc"){
  go_readout <- read.csv(str_glue("~/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/dcHiC_Compartment_Cluster/clustered_dcHIC_PCs_{nature}_k_means_cluster_{cluster}_GO_BP.csv"))
} else if (type == "dc_fuzz"){
  go_readout <- read.csv(str_glue("~/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/Cooltools_compartment_cluster/dcHiC_PC_clustered_cooltools_PCs_sample_wise_{nature}_k_means_cluster_{cluster}_GO_BP.csv"))
} else if (type == "dc_sign"){
    go_readout <- read.csv(str_glue("~/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/Sign_wise_compartment_cluster/dcHiC_PC_Qnm_sign_clustered_cluster_{cluster}_GO_BP.csv"))
} else if (type == "dc_thres"){
  go_readout <- read.csv(str_glue("~/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/Cooltools_compartment_cluster/dcHiC_thres_{thres_decided}_clustered_k_means_cluster_{cluster}_GO_BP.csv"))
}else if (type == "tad"){
  go_readout <- read.csv(str_glue("~/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/Spectral_TAD_calling/Overall_Regeneration_topGO_{cluster}_BP_Ontology.csv"))
} else if (type == "tad_relabelled"){
  go_readout <- read.csv(str_glue("~/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/Spectral_TAD_calling/Overall_Regeneration_topGO_{cluster}_BP_Ontology_relabelled.csv"))
} else if (type == "two_tp_tad"){
  go_readout <- read.csv(str_glue("~/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/Spectral_TAD_calling/TAD_compare/GO_tables/{two_timepoint_selected}_topGO_{cluster}_BP_Ontology.csv"))
}
go_readout <- go_readout %>% mutate(weight_minus_log_10_p = -log10(weight))
go_readout <- go_readout %>% mutate(fraction_annotated_significant = Significant/Annotated)
go_readout <- go_readout %>% arrange(weight_minus_log_10_p)
plt <- ggplot(go_readout) +
  geom_col(aes(Significant,reorder(Term,-weight),fill = weight)) +scale_fill_gradientn(colors = c("red", "purple", "blue")) +xlab("Significant Genes") + ylab("GO Terms") + ggtitle(if (type == "two_tp_tad") {str_glue("GO BP {two_timepoint_selected} {cluster} TAD")} else{str_glue("GO BP {type} cluster {cluster}")}) +labs(fill = "weight p-value")+theme(text = element_text(size = 20))
#scale_fill_distiller(type = "seq",palette = "Spectral", direction = 1) 
plt 
if(type == "cool"){
  ggsave(str_glue("~/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/Cooltools_compartment_cluster/cooltools_{nature}_k_means_cluster_{cluster}_GO_BP_barplot.png"),width = 20, height = 20, units = "cm")
} else if (type == "dc"){
  ggsave(str_glue("~/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/dcHiC_Compartment_Cluster/dcHIC_{nature}_k_means_cluster_{cluster}_GO_BP_barplot.png"),width = 20, height = 20, units = "cm")
} else if (type == "dc_fuzz"){
  ggsave(str_glue("~/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/Cooltools_compartment_cluster/dcHiC_PC_{nature}_k_means_cluster_{cluster}_GO_BP_barplot.png"),width = 20, height = 20, units = "cm")
} else if (type == "dc_sign"){
  ggsave(str_glue("~/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/Sign_wise_compartment_cluster/dcHiC_PC_Qnm_sign_clustered_cluster_{cluster}_GO_BP_barplot.png"),width = 20, height = 20, units = "cm")
} else if (type == "dc_thres"){
  ggsave(str_glue("~/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/Cooltools_compartment_cluster/dcHiC_thres_{thres_decided}_clustered_k_means_cluster_{cluster}_GO_BP_barplot.svg"),width = 40, height = 40, units = "cm")
} else if (type == "tad"){
  ggsave(str_glue("~/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/Spectral_TAD_calling/Regeneration_topGO_{cluster}_GO_BP_barplot.png"),width = 20, height = 20, units = "cm")
} else if (type == "tad_relabelled"){
  ggsave(str_glue("~/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/Spectral_TAD_calling/Relabelled_Regeneration_topGO_{cluster}_GO_BP_barplot.png"),width = 20, height = 20, units = "cm")
} else if (type == "two_tp_tad"){
  ggsave(str_glue("~/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/Spectral_TAD_calling/TAD_compare/GO_tables/{two_timepoint_selected}_topGO_{cluster}_GO_BP_barplot.svg"),width = 40, height = 40, units = "cm")
}
plt <- ggplot(go_readout,aes(x = reorder(Term,-weight),y = Significant)) +
  geom_segment(aes(x = Term,xend = Term,y=0,yend = Significant,color = weight))+scale_color_gradientn(colors = c("red", "purple", "blue")) +ylab("Significant Genes") + xlab("GO Terms") + ggtitle(if (type == "two_tp_tad") {str_glue("GO BP {two_timepoint_selected} {cluster} TAD")} else{str_glue("GO BP {type} cluster {cluster}")}) +labs(color = "weight p-value") +  geom_point(aes(size = fraction_annotated_significant*100),color = "black",alpha = 0.3) + geom_text(aes(label = round(fraction_annotated_significant*100)), color = "yellow", size = 1.00) + 
labs(size = "% of annotated genes significant") +
  coord_flip()#scale_fill_distiller(type = "seq",palette = "Spectral", direction = 1) 
plt 
if(type == "cool"){
  ggsave(str_glue("~/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/Cooltools_compartment_cluster/cooltools_{nature}_k_means_cluster_{cluster}_GO_BP_lollipop_plot.png"),width = 20, height = 20, units = "cm")
} else if (type == "dc"){
  ggsave(str_glue("~/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/dcHiC_Compartment_Cluster/dcHIC_{nature}_k_means_cluster_{cluster}_GO_BP_lollipop_plot.png"),width = 20, height = 20, units = "cm")
} else if (type == "dc_fuzz"){
  ggsave(str_glue("~/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/Cooltools_compartment_cluster/dcHiC_PC_{nature}_k_means_cluster_{cluster}_GO_BP_lollipop_plot.png"),width = 20, height = 20, units = "cm")
} else if (type == "dc_sign"){
  ggsave(str_glue("~/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/Sign_wise_compartment_cluster/dcHiC_PC_Qnm_sign_clustered_cluster_{cluster}_GO_BP_lollipop_plot.png"),width = 20, height = 20, units = "cm")
} else if (type == "dc_thres"){
  ggsave(str_glue("~/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/Cooltools_compartment_cluster/dcHiC_thres_{thres_decided}_clustered_k_means_cluster_{cluster}_GO_BP_lollipop_plot.png"),width = 20, height = 20, units = "cm")
} else if (type == "tad"){
  ggsave(str_glue("~/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/Spectral_TAD_calling/Regeneration_topGO_{cluster}_GO_BP_lollipop_plot.png"),width = 20, height = 20, units = "cm")
} else if (type == "tad_relabelled"){
  ggsave(str_glue("~/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/Spectral_TAD_calling/Relabelled_Regeneration_topGO_{cluster}_GO_BP_lollipop_plot.png"),width = 20, height = 20, units = "cm")
} else if (type == "two_tp_tad"){
  ggsave(str_glue("~/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/Spectral_TAD_calling/TAD_compare/GO_tables/{two_timepoint_selected}_topGO_{cluster}_GO_BP_lollipop_plot.png"),width = 20, height = 20, units = "cm")
}
}

if (type == "tad"){
  TAD_df_overall <- read_feather(str_glue("~/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/Spectral_TAD_calling/Time_TADs_overall.feather"))
  Unique_TAD_cats <- unique(TAD_df_overall$Category)
  map(Unique_TAD_cats,cluster_plotter)
  } else if (type == "tad_relabelled") {
    TAD_df_reread <- read_feather("~/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/Spectral_TAD_calling/Relabelled_Time_TADs_overall.feather")
    Unique_TAD_cats <- unique(TAD_df_reread$Category_relabelled)
    map(Unique_TAD_cats,cluster_plotter)
  } else if (type == "two_tp_tad") {
    Two_tp_df_reread <- read.csv(str_glue("~/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/Spectral_TAD_calling/TAD_compare/{two_timepoint_selected}_50KB_TADs.csv"))
    Unique_TAD_cats <- na.omit(unique(Two_tp_df_reread$Type))
    map(Unique_TAD_cats,cluster_plotter)
  }else{
 map(seq(1,cluster_max,by=1),cluster_plotter)}