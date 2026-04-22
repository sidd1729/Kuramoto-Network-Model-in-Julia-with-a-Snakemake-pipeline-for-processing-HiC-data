Sys.setenv(http_proxy = '172.16.2.252:3128')
Sys.setenv(https_proxy = '172.16.2.252:3128')

library(topGO)
library(dplyr)
library(stringr)
library(rGREAT)
library(GenomicRanges)
library(purrr)
library(enrichR)
library(ggplot2)

listEnrichrSites()

#Setting site as FishEnrichr
setEnrichrSite("FishEnrichr")

listEnrichrDbs()

dbs <- c("GO_Molecular_Function_2018","GO_Cellular_Component_2018","GO_Biological_Process_2018")

chromosome_no <- 25
read_options <- c("control_hpi_12","hpi_12_dpi_2","dpi_2_dpi_4","dpi_4_dpi_7")
chosen_option <- read_options[1]

TAD_compared_2_tp <- readRDS(str_glue("~/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/Spectral_TAD_calling/TAD_compare/{chosen_option}_50KB_TADs.rds"))
All_chrom_TAD_compared_list <- list(TAD_compared_2_tp[[1]]$TAD_Frame)
All_chrom_TAD_compared_list <- rep(All_chrom_TAD_compared_list,25)
TAD_compared_2_tp_1<-TAD_compared_2_tp[[1]]
for(i in seq(1,chromosome_no,by =1)){
TAD_compared_2_tp_chrom <- TAD_compared_2_tp[[i]]$TAD_Frame
TAD_compared_2_tp_chrom <- TAD_compared_2_tp_chrom %>% mutate(Chrom = str_glue("chr{i}"))
All_chrom_TAD_compared_list[[i]] <- TAD_compared_2_tp_chrom
}
Overall_TADs <- do.call(rbind,All_chrom_TAD_compared_list)
write.csv(Overall_TADs,str_glue("~/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/Spectral_TAD_calling/TAD_compare/{chosen_option}_50KB_TADs.csv"))

#Getting background gene set for ontology
two_tp_var_boundary <- Overall_TADs %>% dplyr::select(Boundary,Chrom) %>% mutate(chr = Chrom, start = Boundary, end = Boundary) %>% dplyr::select(chr, start, end)

TAD_range <- makeGRangesFromDataFrame(two_tp_var_boundary)
danRer11_extended<-extendTSSFromOrgDb("danRer11",)
#extendTSSFromDataFrame(time_var_boundary,genome = "danRer11",gene_id_type = "GO:BP")
great_shift <- great(TAD_range, "GO:BP",biomart_dataset = "drerio_gene_ensembl",extended_tss = danRer11_extended )
#job<-great(TAD_range,"GO:BP", biomart_dataset = "drerio_gene_ensembl") #Assume this runs danRer11, not G
gene_assoc_overall = getRegionGeneAssociations(great_shift)
overall_background_gene_set <- unlist(gene_assoc_overall$annotated_genes)
write.csv(overall_background_gene_set,str_glue("~/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/Spectral_TAD_calling/TAD_compare/GO_tables/{chosen_option}_overall_genes.csv"))
Successful_gene_name_go_mapping <- readMappings("~/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/Ontology/zfin_gene_name_go.map")

#Time_comparison_list_ne = list(sparse_mats[[2]],sparse_mats_2[[2]],sparse_mats_3[[2]],sparse_mats_4[[2]],sparse_mats_5[[2]])
#time_var_ne = TimeCompare(Time_comparison_list_ne,resolution = 50000)
#time_var_ne_filt <- time_var_ne$TAD_Bounds

#Getting gene symbols associated with every TAD boundary
Unique_TAD_cats <- na.omit(unique(Overall_TADs$Type))

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
  time_var_ne_filt <- Overall_TADs %>% dplyr::filter(Type == cat)
  time_var_boundary <- time_var_ne_filt %>% dplyr::select(Boundary,Chrom) %>% mutate(chr = Chrom, start = Boundary, end = Boundary) %>% dplyr::select(chr, start, end)
  
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
  write.csv(enrichment_table,file = str_glue("~/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/Spectral_TAD_calling/TAD_compare/GO_tables/{chosen_option}_GREAT_{cat}_BP_Ontology.csv"))
  
  
  #topGO portion
  myInterestingGenes <- unlist(gene_assoc_annotated)
  write.csv(myInterestingGenes,str_glue("~/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/Spectral_TAD_calling/TAD_compare/GO_tables/{chosen_option}_{cat}_genes.csv"))
  #Choosing the background and foreground for the run
  geneNames <- overall_background_gene_set
  #myInterestingGenes <- sample(geneNames, length(geneNames) / 10)
  #myInterestingGenes <- unique(Cleared_c1$ZFIN.ID)
  geneList <- factor(as.integer(geneNames %in% myInterestingGenes))
  if (nlevels(geneList) == 1){
    next
  }
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
  write.csv(allRes,file = str_glue("~/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/Spectral_TAD_calling/TAD_compare/GO_tables/{chosen_option}_topGO_{cat}_BP_Ontology.csv"))
  
  #Fishenrichr portion
  #ensembl = useMart("ensembl",dataset = "drerio_gene_ensembl")
  #gene_name <- getBM(attributes = "external_gene_name" ,filters = 'ensembl_gene_id',values = unlist(gene_assoc$annotated_genes),mart = ensembl)
  
  df <- enrichr(myInterestingGenes, databases = dbs, background = NULL, include_overlap = FALSE, sleepTime = 1)
  write.csv(df$GO_Biological_Process_2018,file = str_glue("~/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/Spectral_TAD_calling/TAD_compare/GO_tables/{chosen_option}_{cat}_fishEnrichR_BP_Ontology.csv"))
  
  #png(str_glue("~/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/Spectral_TAD_calling/TAD_compare/GO_tables/{chosen_option}_{cat}_fishEnrichR_BP_Ontology.png"))
  
  tempPlot1 <- plotEnrich(df$GO_Biological_Process_2018, showTerms = 20, numChar = 60, y = "Count", orderBy = "P.value", xlab = NULL, ylab = NULL, title = str_glue("{chosen_option} {cat} BP"))
  ggsave(tempPlot1,file = str_glue("~/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/Spectral_TAD_calling/TAD_compare/GO_tables/{chosen_option}_{cat}_fishEnrichR_BP_Ontology.png"))
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