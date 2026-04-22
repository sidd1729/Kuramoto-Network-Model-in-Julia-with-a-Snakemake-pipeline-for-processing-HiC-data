Sys.setenv(http_proxy = '172.16.2.252:3128')
Sys.setenv(https_proxy = '172.16.2.252:3128')

#Spectral TAD Clustering
library(SpectralTAD)
library(TADCompare)
library(stringr)
library(dplyr)
library(rGREAT)
library(GenomicRanges)
library(arrow)
#library(enrichplot)
library(biomaRt)
library(enrichR)
library(ggplot2)
library(GOplot)
library(topGO)
library(purrr)
listEnrichrSites()

#Setting site as FishEnrichr
setEnrichrSite("FishEnrichr")

listEnrichrDbs()

dbs <- c("GO_Molecular_Function_2018","GO_Cellular_Component_2018","GO_Biological_Process_2018")
#Databases to query in Fish
#Get the rao contact matrix built into the package
#data("rao_chr20_25_rep")
#head(rao_chr20_25_rep)

#We see that this is a sparse 3-column contact matrix
#Running the algorithm with resolution specified
#results = SpectralTAD(rao_chr20_25_rep, chr = "chr20", resolution = 25000, qual_filter = FALSE, z_clust = FALSE)
#Printing the top 5 TADs
#head(results$Level_1, 5)


#Read in data
#cool_mat = read.table("C:/Users/siddh/Downloads/Rao.GM12878.50kb.txt")
#cool_mat_2 = read.table("C:/Users/siddh/Downloads/50kb_cools/mats/2_dpi_S208_50kb.txt")

#data <- read.table("C:/Users/siddh/Downloads/50kb_cools/mats/2_dpi_S208_50kb.txt", 
                   #encoding = "UTF-16LE", 
                   #skipNul = TRUE, 
                   #header = FALSE,  # Adjust based on your file
                   #sep = "\t")      # Typical for Hi-C matrices

#Convert to sparse 3-column matrix using cooler2sparse from HiCcompare
#sparse_mats = HiCcompare::cooler2sparse(cool_mat)
#Remove empty matrices if necessary
#sparse_mats = sparse_mats$cis[sapply(sparse_mats, nrow) != 0]

#Run SpectralTAD

dataset <- "control_S68"
dataset_2 <- "12_hpi_S207"
dataset_3 <- "2_dpi_S208"
dataset_4 <- "4_dpi_S69"
dataset_5 <- "7_dpi_S1"
#dataset <- dataset_5
sparse_mats = readRDS(str_glue("~/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/temporaries/res_specific_sparse/{dataset}_50kb_sparse.rds"))
#sparse_mats_2 = readRDS("C:/Users/siddh/Downloads/sparse_rds/2_dpi_S208_50kb_sparse.rds")
spec_tads = lapply(names(sparse_mats), function(x) {
  SpectralTAD(sparse_mats[[x]], chr = x,qual_filter = FALSE,levels = 3)
})

saveRDS(spec_tads,str_glue("~/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/temporaries/res_specific_tads/{dataset}_50KB_TADs.rds"))
#spec_tads_2 <- readRDS("C:/Users/siddh/Downloads/spec_tads/2_dpi_S208_50KB_TADs.rds")
#obj_see = spec_tads[[1]][[1]]
#obj_see_2 = spec_tads_2[[1]][[1]]

sparse_mats = readRDS(str_glue("~/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/temporaries/res_specific_sparse/{dataset}_50kb_sparse.rds"))
sparse_mats_2 = readRDS(str_glue("~/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/temporaries/res_specific_sparse/{dataset_2}_50kb_sparse.rds"))
sparse_mats_3 = readRDS(str_glue("~/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/temporaries/res_specific_sparse/{dataset_3}_50kb_sparse.rds"))
sparse_mats_4 = readRDS(str_glue("~/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/temporaries/res_specific_sparse/{dataset_4}_50kb_sparse.rds"))
sparse_mats_5 = readRDS(str_glue("~/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/temporaries/res_specific_sparse/{dataset_5}_50kb_sparse.rds"))

spec_tads = readRDS(str_glue("~/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/temporaries/res_specific_tads/{dataset}_50KB_TADs.rds"))
spec_tads_2 = readRDS(str_glue("~/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/temporaries/res_specific_tads/{dataset_2}_50KB_TADs.rds"))
spec_tads_3 = readRDS(str_glue("~/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/temporaries/res_specific_tads/{dataset_3}_50KB_TADs.rds"))
spec_tads_4 = readRDS(str_glue("~/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/temporaries/res_specific_tads/{dataset_4}_50KB_TADs.rds"))
spec_tads_5 = readRDS(str_glue("~/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/temporaries/res_specific_tads/{dataset_5}_50KB_TADs.rds"))


tad_compare = function(x) {
  chr_bound_1 = bind_rows(spec_tads_4[[x]])
  chr_bound_2 = bind_rows(spec_tads_5[[x]])
  bound_bed <- list(chr_bound_1,chr_bound_2)
  return (TADCompare(sparse_mats_4[[x]],sparse_mats_5[[x]],resolution = 50000, pre_tads = bound_bed))
}
tad_mapped <- map(seq(1,25,by=1),tad_compare)
#saveRDS(tad_mapped,"~/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/Spectral_TAD_calling/TAD_compare/dpi_4_dpi_7_50KB_TADs.rds")

names(sparse_mats)
chr1_bound_1 = bind_rows(spec_tads[[1]])
chr1_bound_2 = bind_rows(spec_tads_2[[1]])
bound_bed <- list(chr1_bound_1,chr1_bound_2)
results = TADCompare(sparse_mats$chr1,sparse_mats_2$chr1,resolution = 50000, pre_tads = bound_bed)

# Visualizing the results
p <- DiffPlot(tad_diff  = results, 
         cont_mat1   = sparse_mats$chr1,
         cont_mat2   = sparse_mats_2$chr1,
         resolution  = 50000,
         start_coord = 500000,
         end_coord   = 700000,
         pre_tad     = bound_bed,
         show_types  = FALSE, 
         point_size  = 5,
         palette     = "RdYlBu",
         rel_heights = c(1, 1))
plot(p)

chr_seq <-seq(1,25,by = 1)
for (i in chr_seq){
  #i=1
Time_comparison_list_1 = list(sparse_mats[[i]],sparse_mats_2[[i]],sparse_mats_3[[i]],sparse_mats_4[[i]],sparse_mats_5[[i]])
time_var <- TimeCompare(Time_comparison_list_1, resolution = 50000)
TAD_df <- time_var$TAD_Bounds
TAD_df <- TAD_df %>% mutate(Chrom = str_glue("chr{i}"))
if (i ==1){
  TAD_df_list <- list(TAD_df)
  rep(TAD_df_list,length(chr_seq))
}
TAD_df_list[[i]] = TAD_df
write_feather(TAD_df,str_glue("~/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/Spectral_TAD_calling/Time_TADs_{i}.feather"))
write.csv(TAD_df,str_glue("~/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/Spectral_TAD_calling/Time_TADs_{i}.csv"))
}

TAD_df_overall <- do.call(rbind,TAD_df_list)
write_feather(TAD_df_overall,str_glue("~/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/Spectral_TAD_calling/Time_TADs_overall.feather"))
write.csv(TAD_df_overall,str_glue("~/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/Spectral_TAD_calling/Time_TADs_overall.csv"))

#Plotting chromosome wise different TAD types
chrom_tad_type_plot <- ggplot(TAD_df_overall, aes(x = Chrom, fill = Category)) + theme(axis.text.x = element_text(size=6)) +
  geom_bar()
chrom_tad_type_plot
ggsave("~/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/Spectral_TAD_calling/Overall_TAD_types_chromwise_barplot.png")

#Plotting TAD types for every chromosome together
TAD_df_overall <- TAD_df_overall %>%
  group_by(Category) %>%
  summarise(counts = n()) #n gives group size in summarise
TAD_df_overall
overall_tad_type_plot <- ggplot(TAD_df_overall, aes(x = Category, y = counts,fill = Category)) + theme(axis.text.x = element_text(size=6)) +
  geom_bar(stat = "identity")#,position = "stack")
overall_tad_type_plot
ggsave("~/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/Spectral_TAD_calling/Overall_TAD_types_aggregated_barplot.png")



TAD_df_overall <- do.call(rbind,TAD_df_list)
#Getting background gene set for ontology
time_var_boundary <- TAD_df_overall %>% dplyr::select(Coordinate,Chrom) %>% mutate(chr = Chrom, start = Coordinate, end = Coordinate) %>% dplyr::select(chr, start, end)

TAD_range <- makeGRangesFromDataFrame(time_var_boundary)
danRer11_extended<-extendTSSFromOrgDb("danRer11",)
#extendTSSFromDataFrame(time_var_boundary,genome = "danRer11",gene_id_type = "GO:BP")
great_shift <- great(TAD_range, "GO:BP",biomart_dataset = "drerio_gene_ensembl",extended_tss = danRer11_extended )
#job<-great(TAD_range,"GO:BP", biomart_dataset = "drerio_gene_ensembl") #Assume this runs danRer11, not G
gene_assoc_overall = getRegionGeneAssociations(great_shift)
overall_background_gene_set <- unlist(gene_assoc_overall$annotated_genes)
write.csv(overall_background_gene_set,str_glue("~/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/Spectral_TAD_calling/Overall_Regeneration_Time_TAD_genes.csv"))
Successful_gene_name_go_mapping <- readMappings("~/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/Ontology/zfin_gene_name_go.map")

#Time_comparison_list_ne = list(sparse_mats[[2]],sparse_mats_2[[2]],sparse_mats_3[[2]],sparse_mats_4[[2]],sparse_mats_5[[2]])
#time_var_ne = TimeCompare(Time_comparison_list_ne,resolution = 50000)
#time_var_ne_filt <- time_var_ne$TAD_Bounds

#Getting gene symbols associated with every TAD boundary
Unique_TAD_cats <- unique(TAD_df_overall$Category)

plot_list_BP = list()
#plot_list_CC = list()
#plot_list_MF = list()

##Function to find significant GO terms taken from https://stackoverflow.com/questions/79504284/how-to-extract-the-significant-genes-from-a-go-analysis-with-topgo-in-r (Answer by Gwang-Jin Kim)
# Which are the significant GO terms?
significantGO <- function(result_test, alpha = 0.01) {
  if (!inherits(result_test, "topGOresult")) {
    stop("result_test must be a topGOresult object.")
  }
  
  significant_terms <- names(score(result_test)[score(result_test) < alpha])
  
  if (length(significant_terms) == 0) {
    message("No significant GO terms found at alpha = ", alpha)
  }
  
  return(significant_terms)
}

for (cat in Unique_TAD_cats){
#cat <- Unique_TAD_cats[1]
time_var_ne_filt <- TAD_df_overall %>% dplyr::filter(Category == cat)
time_var_boundary <- time_var_ne_filt %>% dplyr::select(Coordinate,Chrom) %>% mutate(chr = Chrom, start = Coordinate, end = Coordinate) %>% dplyr::select(chr, start, end)

#con_tads <- ConsensusTADs(time_mats, resolution = 50000)
#head(con_tads$Consensus)

#TAD_Frame <- results$TAD_Frame
#TAD_Frame <- TAD_Frame %>% dplyr::filter((Type == "Strength Change") & 
                                           #(Enriched_In == "Matrix 2"))

# Assign a chromosome and convert to a bed format
#TAD_Frame <- TAD_Frame %>% dplyr::select(Boundary) %>% mutate(chr = "chr1", 
                                                              #start = Boundary, end = Boundary) %>% dplyr::select(chr, start, end)

# Set up rGREAT job with default parameters
#TAD_range <- makeGRangesFromDataFrame(TAD_Frame)
TAD_range <- makeGRangesFromDataFrame(time_var_boundary)


#great_shift <- great(TAD_Frame, "GO:BP", "TxDb.Drerio.UCSC.danRer11.refGene")
job<-great(TAD_range,"GO:BP", biomart_dataset = "drerio_gene_ensembl",extended_tss = danRer11_extended) #Assume this runs danRer11, not G
#job <- great_shift
#job = great(TAD_Frame, 
            #geneSets = "GO:BP", 
            #species = "danRer11", 
            #mode = "basalPlusExt")

# Get results
enrichment_table = getEnrichmentTable(job)
gene_assoc = getRegionGeneAssociations(job)
gene_assoc_annotated <- gene_assoc$annotated_genes
#enrichment_table <- enrichment_table %>% mutate(corres_gene = I(unlist(gene_assoc_annotated)))

#plotRegionGeneAssociations(job)
#plotEnrich(enrichment_table,showTerms = 20, numChar = 60, y = "observed_region_hits")
write.csv(enrichment_table,file = str_glue("~/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/Spectral_TAD_calling/Overall_Regeneration_GREAT_{cat}_BP_Ontology.csv"))


#topGO portion
myInterestingGenes <- unlist(gene_assoc_annotated)
write.csv(myInterestingGenes,str_glue("~/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/Spectral_TAD_calling/Overall_Regeneration_{cat}_genes.csv"))
#Choosing the background and foreground for the run
geneNames <- overall_background_gene_set
#myInterestingGenes <- sample(geneNames, length(geneNames) / 10)
#myInterestingGenes <- unique(Cleared_c1$ZFIN.ID)
geneList <- factor(as.integer(geneNames %in% myInterestingGenes))
names(geneList) <- geneNames
#str(geneList)
#Maybe creating a topGO object?
 GOdata <- new("topGOdata", ontology = "BP", allGenes = geneList,
                annot = annFUN.gene2GO, gene2GO = Successful_gene_name_go_mapping)
#relevant_genes <- genesInTerm(GOdata)
#strange <- str(relevant_genes)
test.stat <- new("classicScore", 
                  testStatistic = GOKSTest, 
                  name = "KS tests")
  
  resultKS <- getSigGroups(GOdata, test.stat)
  test.stat <- new("weightCount", testStatistic = GOFisherTest, 
                   name = "Fisher test", sigRatio = "ratio")
  resultWeight <- getSigGroups(GOdata, test.stat)
  resultFis <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
  allRes <- GenTable(GOdata, classic = resultFis, KS = resultKS, 
                     weight = resultWeight, orderBy = "weight", 
                     ranksOf = "classic", topNodes = 20,numChar = 200)
  
  ##Getting genes with significant association with GO terms using the aforementioned stackoverflow answer
  significant_GO_terms <- significantGO(resultWeight, alpha = max(as.numeric(max(allRes$weight)) + 0.001,0.05)) 
  
  # Which genes are associated with these GO terms?
  significant_genes_list <- genesInTerm(GOdata,whichGO = significant_GO_terms)
  
  # Format nicely to Gene - GO Term pairs
  significant_genes_df <- stack(significant_genes_list)
  colnames(significant_genes_df) <- c("Gene", "GO.ID")
  
  # Keep only the significant genes
  significant_genes_df <- significant_genes_df %>% filter(Gene %in% names(geneList[geneList == 1]))
  
  #significantly_clean <- significant_genes_df %>%  distinct(GO.ID,Gene, .keep_all = TRUE)
  relevant_genes <-  split(significant_genes_df$Gene,significant_genes_df$GO.ID)
  
  #toString(relevant_genes[["GO:0015671"]])
  GO_IDs_mapped<-map(allRes$GO.ID, function(x) relevant_genes[[str_glue("{x}")]])
  allRes <- allRes %>% mutate(gene_names = I(GO_IDs_mapped) )
  allRes$gene_names   <- unlist(map(allRes$gene_names, function(x) toString(x)))
  allRes$gene_names  <- unlist(map(allRes$gene_names, function(x) sapply(strsplit(x, ","), paste, collapse=" ")))
  write.csv(allRes,file = str_glue("~/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/Spectral_TAD_calling/Overall_Regeneration_topGO_{cat}_BP_Ontology.csv"))
  
#Fishenrichr portion
#ensembl = useMart("ensembl",dataset = "drerio_gene_ensembl")#,host = "www.ensembl.org")
#gene_name <- getBM(attributes = "external_gene_name" ,filters = 'ensembl_gene_id',values = unlist(gene_assoc$annotated_genes),mart = ensembl)
#gene_name <- getBM(attributes = "external_gene_name" ,filters = 'ensembl_gene_id',values = unlist(gene_assoc$annotated_genes),mart = ensembl)

#df <- enrichr(gene_name$external_gene_name, databases = dbs, background = NULL, include_overlap = FALSE, sleepTime = 1)
df <- enrichr(myInterestingGenes, databases = dbs, background = NULL, include_overlap = FALSE, sleepTime = 1)
write.csv(df$GO_Biological_Process_2018,file = str_glue("~/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/Spectral_TAD_calling/Overall_Regeneration_{cat}_fishEnrichR_BP_Ontology.csv"))

#png(str_glue("~/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/Spectral_TAD_calling/Overall_Regeneration_{cat}_fishEnrichR_BP_Ontology.png"))

tempPlot1 <- plotEnrich(df$GO_Biological_Process_2018, showTerms = 20, numChar = 60, y = "Count", orderBy = "P.value", xlab = NULL, ylab = NULL, title = str_glue("{cat} BP"))
ggsave(tempPlot1,file = str_glue("~/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/Spectral_TAD_calling/Overall_Regeneration_{cat}_fishEnrichR_BP_Ontology.png"))
#dev.off()


#Fishenrichr GO BP and other parts
#write.csv(df$GO_Molecular_Function_2018,file = str_glue("~/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/Spectral_TAD_calling/Overall_Regeneration_{cat}_fishEnrichR_MF_Ontology.csv"))

#png(str_glue("~/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/Spectral_TAD_calling/Overall_Regeneration_{cat}_fishEnrichR_MF_Ontology.png"))

#tempPlot2 <- plotEnrich(df$GO_Molecular_Function_2018, showTerms = 20, numChar = 60, y = "Count", orderBy = "P.value", xlab = NULL, ylab = NULL, title = str_glue("{cat}_MF"))

#ggsave(tempPlot2,file = str_glue("~/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/Spectral_TAD_calling/Overall_Regeneration_{cat}_fishEnrichR_MF_Ontology.png"))
#dev.off()

#write.csv(df$GO_Cellular_Component_2018,file = str_glue("~/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/Spectral_TAD_calling/Overall_Regeneration_{cat}_fishEnrichR_CC_Ontology.csv"))

#png(str_glue("~/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/Spectral_TAD_calling/Overall_Regeneration_{cat}_fishEnrichR_CC_Ontology.png"))

#tempPlot3 <- plotEnrich(df$GO_Cellular_Component_2018, showTerms = 20, numChar = 60, y = "Count", orderBy = "P.value", xlab = NULL, ylab = NULL, title = str_glue("{cat}_CC"))

#ggsave(tempPlot3,file = str_glue("~/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/Spectral_TAD_calling/Overall_Regeneration_{cat}_fishEnrichR_CC_Ontology.png"))

#dev.off()
}
# Subset to only include vital information
#enrichment_table <- bind_rows(enrichment_table, .id = "source") %>% 
  #dplyr::select(Ontology = source, Description = name, 
                #`P-value` = Hyper_Raw_PValue)

# Print head organizaed by p-values
#head(enrichment_table %>% dplyr::arrange(`P-value`))

#test_table<-read_feather("~/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/Spectral_TAD_calling/filtered_diff_compartments.feather")
#test_table_chr_1 <- test_table %>% filter(chr == "chr1")
#test_grange_1 <- makeGRangesFromDataFrame(test_table_chr_1)
#job_diff<-great(test_grange_1,"GO:BP", biomart_dataset = "drerio_gene_ensembl") #Assume this runs danRer11, not G
#job = great(TAD_Frame, 
#geneSets = "GO:BP", 
#species = "danRer11", 
#mode = "basalPlusExt")

# Get results
#enrichment_table_diff_1 = getEnrichmentTable(job_diff)
#write.csv(enrichment_table_diff_1,file = "~/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/Spectral_TAD_calling/Chr1_Regeneration_Differential_Compartment_BP_Ontology.csv")
