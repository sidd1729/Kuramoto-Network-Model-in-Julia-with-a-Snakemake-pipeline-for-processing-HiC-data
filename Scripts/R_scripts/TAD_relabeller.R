Sys.setenv(http_proxy = '172.16.2.252:3128')
Sys.setenv(https_proxy = '172.16.2.252:3128')

#Reading in a timeCompare output first
library(dplyr)
library(rGREAT)
library(topGO)
library(arrow)
library(GenomicRanges)
library(stringr)
library(ggplot2)
library(purrr)
library(enrichR)
listEnrichrSites()

#Setting site as FishEnrichr
setEnrichrSite("FishEnrichr")

listEnrichrDbs()

dbs <- c("GO_Molecular_Function_2018","GO_Cellular_Component_2018","GO_Biological_Process_2018")

TAD_Compare_reread <- read.csv("~/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/Spectral_TAD_calling/Time_TADs_overall.csv")
TAD_Compare_reread_slice <- TAD_Compare_reread[c(3,4,5,6,7)]
Three_check <- apply(TAD_Compare_reread_slice,2,function(x) x>=3) #This seems to be correct margin for rowwise operations on df
###
#Timepoint TAD counting
tp_1_boundaries <- sum(Three_check[,1])
tp_2_boundaries <- sum(Three_check[,2])
tp_3_boundaries <- sum(Three_check[,3])
tp_4_boundaries <- sum(Three_check[,4])
tp_5_boundaries <- sum(Three_check[,5])

variance_check <- function(x){
  result = case_when(
    (x[1] & x[2]) == TRUE ~ "Invariant",
     (!x[1]) & (x[2]) == TRUE | (x[1]) & (!x[2]) == TRUE ~ "Variant"
  )
  return(result)
}
variance_check_1_2 <- apply(Three_check[,c(1,2)],1,variance_check)
variants_1_2 <- sum(variance_check_1_2 == "Variant",na.rm = TRUE)
invariants_1_2 <- sum(variance_check_1_2 == "Invariant",na.rm = TRUE)
variance_check_2_3 <- apply(Three_check[,c(2,3)],1,variance_check)
variants_2_3 <- sum(variance_check_2_3 == "Variant",na.rm = TRUE)
invariants_2_3 <- sum(variance_check_2_3 == "Invariant",na.rm = TRUE)
variance_check_3_4 <- apply(Three_check[,c(3,4)],1,variance_check)
variants_3_4 <- sum(variance_check_3_4 == "Variant",na.rm = TRUE)
invariants_3_4 <- sum(variance_check_3_4 == "Invariant",na.rm = TRUE)
variance_check_4_5 <- apply(Three_check[,c(4,5)],1,variance_check)
variants_4_5 <- sum(variance_check_4_5 == "Variant",na.rm = TRUE)
invariants_4_5 <- sum(variance_check_4_5 == "Invariant",na.rm = TRUE)
Complete_invariants <- (variance_check_1_2 == "Invariant") & (variance_check_2_3 == "Invariant")  & (variance_check_3_4 == "Invariant") & (variance_check_4_5 == "Invariant")
Complete_invariant <- sum(Complete_invariants)
tp_changes_vec <-c("Timepoint_1 total",str_glue("{tp_1_boundaries}"),"Timepoint_2 total",str_glue("{tp_2_boundaries}"),"Timepoint_3_total",str_glue("{tp_3_boundaries}"),"Timepoint_4_total",str_glue("{tp_4_boundaries}"),"Timepoint_5_total",str_glue("{tp_5_boundaries}"),"Complete invariants",str_glue("{Complete_invariant}"),
                                                                                                                                                                  "variants_1_2",str_glue("{variants_1_2}"),"variants_2_3",str_glue("{variants_2_3}"),"variants_3_4",str_glue("{variants_3_4}"),"variants_4_5",str_glue("{variants_4_5}"),
                                                                                                                                                                  "invariants_1_2",str_glue("{invariants_1_2}"),"invariants_2_3",str_glue("{invariants_2_3}"),"invariants_3_4",str_glue("{invariants_3_4}"),"invariants_4_5",str_glue("{invariants_4_5}"))
write.table(data.frame(tp_changes = tp_changes_vec), 
            file = "~/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/Spectral_TAD_calling/tp_specific_tad_changes.txt", 
            sep = "\t",        
            row.names = FALSE, 
            quote = FALSE)


##
inefficient_recategorization <- function(x){
  result = case_when(
    all(c(x[1],x[2],x[3],x[4],x[5]) == c(TRUE,TRUE,TRUE,TRUE,TRUE)) == TRUE | all(c(x[1],x[2],x[3],x[4],x[5]) == c(TRUE,FALSE,TRUE,TRUE,TRUE)) == TRUE | 
      all(c(x[1],x[2],x[3],x[4],x[5]) == c(TRUE,TRUE,FALSE,TRUE,TRUE)) == TRUE | all(c(x[1],x[2],x[3],x[4],x[5]) == c(TRUE,TRUE,TRUE,FALSE,TRUE)) == TRUE ~ "Highly Common TAD",
    all(c(x[1],x[2],x[3],x[4],x[5]) == c(FALSE,TRUE,TRUE,TRUE,TRUE)) == TRUE ~ "Early Appearing TAD",
    all(c(x[1],x[2],x[3],x[4],x[5]) == c(TRUE,FALSE,FALSE,FALSE,FALSE)) == TRUE ~ "Early Disappearing TAD",
    all(c(x[1],x[2],x[3],x[4],x[5]) == c(FALSE,FALSE,TRUE,TRUE,TRUE)) == TRUE ~ "Late Appearing TAD",
    all(c(x[1],x[2],x[3],x[4],x[5]) == c(TRUE,TRUE,FALSE,FALSE,FALSE)) == TRUE | all(c(x[1],x[2],x[3],x[4],x[5]) == c(TRUE,TRUE,TRUE,FALSE,FALSE)) == TRUE ~ "Late Disappearing TAD",
    .default = "Dynamic TAD"
  )
  return(result)
}
New_label <- apply(Three_check,1,inefficient_recategorization) #This seems to be correct margin for rowwise operations on matrix
TAD_Compare_reread <- TAD_Compare_reread %>% mutate(Category_relabelled = New_label)
#I guess from underlying difference in internal representation as stated in the below stackoverflow answer 
#(https://stackoverflow.com/questions/34509578/difference-between-data-frame-and-matrix-indexing)

write.csv(TAD_Compare_reread,"~/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/Spectral_TAD_calling/Relabelled_Time_TADs_overall.csv")
write_feather(TAD_Compare_reread,"~/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/Spectral_TAD_calling/Relabelled_Time_TADs_overall.feather")

####
#Plotting chromosome wise different TAD types
chrom_tad_type_plot <- ggplot(TAD_Compare_reread, aes(x = Chrom, fill = Category_relabelled)) + theme(axis.text.x = element_text(size=6)) +
  geom_bar()
chrom_tad_type_plot
ggsave("~/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/Spectral_TAD_calling/Relabelled_Overall_TAD_types_chromwise_barplot.png")

#Plotting TAD types for every chromosome together
TAD_Compare_reread_copy <- TAD_Compare_reread %>%
  group_by(Category_relabelled) %>%
  summarise(counts = n()) #n gives group size in summarise
TAD_Compare_reread
TAD_Compare_reread_copy
overall_tad_type_plot <- ggplot(TAD_Compare_reread_copy, aes(x = Category_relabelled, y = counts,fill = Category_relabelled)) + theme(axis.text.x = element_text(size=6)) +
  geom_bar(stat = "identity")#,position = "stack")
overall_tad_type_plot
ggsave("~/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/Spectral_TAD_calling/Relabelled_Overall_TAD_types_aggregated_barplot.png")


####
#Getting background gene set for ontology
time_var_boundary <- TAD_Compare_reread %>% dplyr::select(Coordinate,Chrom) %>% mutate(chr = Chrom, start = Coordinate, end = Coordinate) %>% dplyr::select(chr, start, end)

TAD_range <- makeGRangesFromDataFrame(time_var_boundary)
danRer11_extended<-extendTSSFromOrgDb("danRer11",)
#extendTSSFromDataFrame(time_var_boundary,genome = "danRer11",gene_id_type = "GO:BP")
great_shift <- great(TAD_range, "GO:BP",biomart_dataset = "drerio_gene_ensembl",extended_tss = danRer11_extended )
#job<-great(TAD_range,"GO:BP", biomart_dataset = "drerio_gene_ensembl") #Assume this runs danRer11, not G
gene_assoc_overall = getRegionGeneAssociations(great_shift)
overall_background_gene_set <- unlist(gene_assoc_overall$annotated_genes)
write.csv(overall_background_gene_set,str_glue("~/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/Spectral_TAD_calling/Overall_Regeneration_Time_TAD_genes_relabelled.csv"))
Successful_gene_name_go_mapping <- readMappings("~/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/Ontology/zfin_gene_name_go.map")

#Time_comparison_list_ne = list(sparse_mats[[2]],sparse_mats_2[[2]],sparse_mats_3[[2]],sparse_mats_4[[2]],sparse_mats_5[[2]])
#time_var_ne = TimeCompare(Time_comparison_list_ne,resolution = 50000)
#time_var_ne_filt <- time_var_ne$TAD_Bounds

#Getting gene symbols associated with every TAD boundary
Unique_TAD_cats <- unique(TAD_Compare_reread$Category_relabelled)

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
  time_var_ne_filt <- TAD_Compare_reread %>% dplyr::filter(Category_relabelled == cat)
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
  write.csv(enrichment_table,file = str_glue("~/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/Spectral_TAD_calling/Overall_Regeneration_GREAT_{cat}_BP_Ontology_relabelled.csv"))
  
  
  #topGO portion
  myInterestingGenes <- unlist(gene_assoc_annotated)
  write.csv(myInterestingGenes,str_glue("~/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/Spectral_TAD_calling/Overall_Regeneration_{cat}_genes_relabelled.csv"))
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
  write.csv(allRes,file = str_glue("~/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/Spectral_TAD_calling/Overall_Regeneration_topGO_{cat}_BP_Ontology_relabelled.csv"))
  
  #Fishenrichr portion
  #ensembl = useMart("ensembl",dataset = "drerio_gene_ensembl")
  #gene_name <- getBM(attributes = "external_gene_name" ,filters = 'ensembl_gene_id',values = unlist(gene_assoc$annotated_genes),mart = ensembl)
  
  df <- enrichr(myInterestingGenes, databases = dbs, background = NULL, include_overlap = FALSE, sleepTime = 1)
  write.csv(df$GO_Biological_Process_2018,file = str_glue("~/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/Spectral_TAD_calling/Overall_Regeneration_{cat}_fishEnrichR_BP_Ontology_relabelled.csv"))
  
  #png(str_glue("~/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/Spectral_TAD_calling/Relabelled_Overall_Regeneration_{cat}_fishEnrichR_BP_Ontology.png"))
  
  tempPlot1 <- plotEnrich(df$GO_Biological_Process_2018, showTerms = 20, numChar = 60, y = "Count", orderBy = "P.value", xlab = NULL, ylab = NULL, title = str_glue("{cat} BP"))
  ggsave(tempPlot1,file = str_glue("~/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/Spectral_TAD_calling/Relabelled_Overall_Regeneration_{cat}_fishEnrichR_BP_Ontology.png"))
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