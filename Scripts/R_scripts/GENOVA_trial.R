#Trial of GENOVA
library(GENOVA)
tested_mat<-load_contacts("~/Siddharth/Program_Directory/Data/Actual_Data/Raw_cool/40kb/2_dpi_S166.cool",balancing = FALSE)
hic_matrixplot(tested_mat,chrom = "chr1:30,000,000-50,000,000")
subset_mat <- subset(tested_mat,'1')
insulated_scores <- insulation_score(tested_mat,window = 5)#,norm_to = c("chromosome"),norm_fun = log2overmedian)
x <- select_subset(tested_mat, "1", 1, 1589)
image(x)