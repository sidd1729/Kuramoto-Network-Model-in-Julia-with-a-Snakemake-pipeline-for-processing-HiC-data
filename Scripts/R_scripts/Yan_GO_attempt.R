Sys.setenv(http_proxy = '172.16.2.252:3128')
Sys.setenv(https_proxy = '172.16.2.252:3128')

library(biomaRt)
library(arrow)
library(furrr)
library(purrr)
library(topGO)
library(GenomicRanges)
library(stringr)
library(dplyr)
#library(clusterProfiler)
library(ggplot2)
nature <- "raw"
thres_decided = 0.2592
#PC_type <- "cool"
annotation_df <- read.csv("~/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/Ontology/zfin.gaf",sep = '\t',comment.char = "!",header = FALSE)
#query_df <- read_feather(str_glue("~/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/dcHiC_Compartment_Cluster/clustered_dcHIC_PCs_sample_wise_{nature}_k_means_df.feather"))
##query_df <- read_feather(str_glue("~/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/Cooltools_compartment_cluster/dcHiC_PC_clustered_cooltools_PCs_sample_wise_{nature}_k_means_df.feather"))
query_df <- read_feather(str_glue("~/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/Cooltools_compartment_cluster/dcHiC_thres_{thres_decided}_clustered_compartments_yan.feather"))
#query_df <- read_feather(str_glue("~/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/Cooltools_compartment_cluster/clustered_cooltools_PCs_sample_wise_{nature}_k_means_df.feather"))
#query_df$chrom <- unlist(map(query_df$chrom,function(x) str_glue("chr{x}"))) # Have to make it chr x for intersection to be found properly
all_zfin_ids_corres_regions <- read.csv("~/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/Ontology/mart_export.txt",sep = '\t')
all_zfin_ids_corres_regions$Chromosome.scaffold.name <- unlist(map(all_zfin_ids_corres_regions$Chromosome.scaffold.name,function(x) str_glue("chr{x}")))
#Processing annotation files
names(annotation_df) <- c("Database_Designation","Marker_ID","Gene_Symbol","Qualifiers","GO_Term_ID",
                    "Reference_ID","GO_Evidence_Code","Inferred_From","Ontology","Marker_Name","Marker_Synonyms",
"Marker_Type","Taxon","Modification_Date","Assigned_By","Annotation_Extension","Gene_Product_Form_ID")

annotation_df_bp = annotation_df[annotation_df$Ontology == "P",]
Unique_Z_IDs <- unique(annotation_df$Marker_ID)
annotation_filtered_list_len <- length(Unique_Z_IDs)

for (i in seq(1,annotation_filtered_list_len,by=1)){
annotation_df_filtrd = annotation_df_bp[annotation_df_bp$Marker_ID == Unique_Z_IDs[i],]$GO_Term_ID
annotation_df_gene_name_filtered = annotation_df_bp[annotation_df_bp$Marker_ID == Unique_Z_IDs[i],]$Gene_Symbol
if (i == 1){
  annotation_filtered_list <- list(toString(annotation_df_filtrd))
  annotation_gene_filtered_list <- list(annotation_df_gene_name_filtered)
  annotation_filtered_list <- rep(annotation_filtered_list,annotation_filtered_list_len)
  annotation_gene_filtered_list <- rep(annotation_gene_filtered_list,annotation_filtered_list_len)}
else{
  annotation_filtered_list[[i]] <- toString(annotation_df_filtrd)
  annotation_gene_filtered_list[[i]] <- annotation_df_gene_name_filtered
}
}
#Creating a gene name to GO annotation
Unique_Gene_IDs <- unique(annotation_df$Gene_Symbol)
annotation_filtered_gene_list_len <- length(Unique_Gene_IDs)
for (i in seq(1,annotation_filtered_gene_list_len,by=1)){
  annotation_df_gene_GO_filtered = annotation_df_bp[annotation_df_bp$Gene_Symbol == Unique_Gene_IDs[i],]$GO_Term_ID
  if (i == 1){
    annotation_gene_GO_filtered_list <- list(annotation_df_gene_GO_filtered)
    annotation_gene_GO_filtered_list <- rep(annotation_gene_GO_filtered_list,annotation_filtered_gene_list_len)}
    #annotation_gene_filtered_list <- rep(annotation_gene_filtered_list,annotation_filtered_list_len)}
  else{
    annotation_gene_GO_filtered_list[[i]] <- annotation_df_gene_GO_filtered
  }
}
#Getting gene symbol to GO map to use after applying GREAT and seeing nearby regions/associated genes with TAD boundary
#adjacent regions
gene_name_go_mapping_df <- data.frame(gene_symbol = Unique_Gene_IDs,go = I(annotation_gene_GO_filtered_list))

#Getting zfin ID to gene Symbol mapping using the readMapping function (not necessary for ontology)
Unique_annotated_gene_names<-map(annotation_gene_filtered_list,unique)
gene_name_mapping_df <- data.frame(z_id = Unique_Z_IDs,go = I(Unique_annotated_gene_names))

#zfin ID and GO maps
mapping_df <- data.frame(z_id = Unique_Z_IDs,go = I(annotation_filtered_list))
entire_zgene_ids <- data.frame(z_id = Unique_Z_IDs)
#mapping_df$go <- map(mapping_df$go,toString)
write.table(gene_name_go_mapping_df,file = "~/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/Ontology/zfin_gene_name_go.map",sep = '\t',col.names = FALSE,row.names = FALSE,quote = FALSE)
write.table(gene_name_mapping_df,file = "~/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/Ontology/zfin_gene_name.map",sep = '\t',col.names = FALSE,row.names = FALSE,quote = FALSE)
write.table(mapping_df,file = "~/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/Ontology/zfin_GO.map",sep = '\t',col.names = FALSE,row.names = FALSE,quote = FALSE)
write.table(entire_zgene_ids,file = "~/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/Ontology/all_zfin_gene_IDs.tsv",sep = '\t')


#Annotation successfully processed as a map between Zfin IDs and GOs
Successful_map <- readMappings("~/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/Ontology/zfin_GO.map")
Successful_gene_map <- readMappings("~/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/Ontology/zfin_gene_name.map")
complete_gene_universe <- entire_zgene_ids$z_id

#GRanges from dataframes for all genes mapped to their chromosomal location
#As well as all compartments in their clusters in the query df
#chr for dchic chrom for cooltools
All_genes_GRange<-makeGRangesFromDataFrame(all_zfin_ids_corres_regions,keep.extra.columns = TRUE,seqnames.field = "Chromosome.scaffold.name",start.field = "Gene.start..bp.",end.field = "Gene.end..bp.")
All_cluster_GRange <- makeGRangesFromDataFrame(query_df,keep.extra.columns = TRUE,seqnames.field = "chr",start.field = "start", end.field = "end") 
#Cluster1_GRange <- makeGRangesFromDataFrame(query_df[query_df$cluster == 2,],keep.extra.columns = TRUE,seqnames.field = "chr",start.field = "start", end.field = "end")
#,query_df = query_df,All_genes_GRange = All_genes_GRange
Cluster_id_mapper <- function(i,queried_df = query_df,genes_grange = All_genes_GRange){
  Cluster_GRange <- makeGRangesFromDataFrame(queried_df[queried_df$cluster == i,],keep.extra.columns = TRUE,seqnames.field = "chr",start.field = "start", end.field = "end")
  cluster_overlap_obj <- findOverlaps(Cluster_GRange,genes_grange)
  Cleared <- All_genes_GRange[subjectHits(cluster_overlap_obj)]
  return(unique(Cleared$ZFIN.ID))
}

plan(multisession,workers = 8)
All_cluster_ids_list <- future_map(seq(1,max(query_df$cluster),by = 1), Cluster_id_mapper)

#All_cluster_ids_list_nonscaled <- ALL_cluster_ids_list

#Finding overlap of query with subject in subject and using subjectHits to index the hits/matching rows in the subject
overlap_obj<-findOverlaps(All_cluster_GRange,All_genes_GRange)
#c1_overlap_obj <- findOverlaps(Cluster1_GRange,All_genes_GRange)

Cleared <- All_genes_GRange[subjectHits(overlap_obj)]
#Cleared_c1 <- All_genes_GRange[subjectHits(c1_overlap_obj)]

#Cleared_query <- All_cluster_GRange[queryHits(overlap_obj)]
gene_ids_clusters <- unique(Cleared$ZFIN.ID)
viable_ids  <- intersect(gene_ids_clusters,complete_gene_universe)
#write.csv(gene_ids_clusters,str_glue("~/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/Cooltools_compartment_cluster/dcHiC_PC_clustered_cooltools_PCs_sample_wise_{nature}_k_means_cluster_all_genes.csv"))
write.csv(gene_ids_clusters,str_glue("~/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/Cooltools_compartment_cluster/dcHiC_thres_{thres_decided}_clustered_all_genes.csv"))
#all_gene_attempt <- factor(rep(1,length(viable_ids)))
#all_gene_attempt[2] <- 2
#all_gene_attempt[3] <- 2
#names(all_gene_attempt) <- viable_ids

#Choosing the background and foreground for the run
geneNames <- gene_ids_clusters

#myInterestingGenes <- sample(geneNames, length(geneNames) / 10)
#myInterestingGenes <- unique(Cleared_c1$ZFIN.ID)

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

for (i in seq(1,max(query_df$cluster),by = 1)){
#i = 1
myInterestingGenes <- All_cluster_ids_list[[i]]
write.csv(myInterestingGenes,str_glue("~/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/Cooltools_compartment_cluster/dcHiC_thres_{thres_decided}_clustered_k_means_cluster_{i}_genes.csv"))
##write.csv(myInterestingGenes,str_glue("~/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/Cooltools_compartment_cluster/dcHiC_PC_clustered_cooltools_PCs_sample_wise_{nature}_k_means_cluster_{i}_genes.csv"))         
#write.csv(myInterestingGenes,str_glue("~/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/dcHiC_Compartment_Cluster/clustered_dcHIC_PCs_{nature}_k_means_cluster_{i}_genes.csv"))
#write.csv(myInterestingGenes,str_glue("~/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/Cooltools_compartment_cluster/clustered_cooltools_PCs_sample_wise_{nature}_k_means_cluster_{i}_genes.csv"))
geneList <- factor(as.integer(geneNames %in% myInterestingGenes))
names(geneList) <- geneNames
#str(geneList)
#Maybe creating a topGO object?
GOdata <- new("topGOdata", ontology = "BP", allGenes = geneList,
              annot = annFUN.gene2GO, gene2GO = Successful_map)
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

# This is the result
#print(significant_genes_df)
#relevant_genes <- significant_genes_df$

#toString(relevant_genes[["GO:0015671"]])
GO_IDs_mapped<-map(allRes$GO.ID, function(x) relevant_genes[[str_glue("{x}")]])
#Below function takes the list  of zfin IDs corresponding to go IDs, first we take the list and essentially see all
#zfin ids corresponding to a GO ID using 1st map/outer map, inner map takes these zfin ids and maps the gene symbol using
#zfin id to gene symbol map which was previously created, unlisting the lists created by inner map (otherwise each zfin id becomes a list item which isn't what we want, 
#rather we would like all gene symbols corresponding to a GO symbol to become a list item with no further levels, which is what this achieves)
gene_names_mapped <- map(GO_IDs_mapped, function(x) unlist(map(x, function(x) Successful_gene_map[[str_glue("{x}")]])))
#GO_IDs_gene_mapped <- map(GO_IDs_mapped, function(x) x)
#relevant_genes <- genesInTerm(GOdata,allRes$GO.ID)
#allResT <- allRes %>% mutate(test = as.character(GO.ID))
#allResT <- allResT %>% mutate(trial = relevant_genes$GO.ID)
#allRes <- allRes %>% mutate(genes = toString(relevant_genes[as.character(GO.ID)]))
allRes <- allRes %>% mutate(zfin_gene_ids = I(GO_IDs_mapped))
allRes$zfin_gene_ids   <- unlist(map(allRes$zfin_gene_ids, function(x) toString(x)))
allRes$zfin_gene_ids  <- unlist(map(allRes$zfin_gene_ids, function(x) sapply(strsplit(x, ","), paste, collapse=" ")))
allRes <- allRes %>% mutate(gene_names = I(gene_names_mapped) )
allRes$gene_names   <- unlist(map(allRes$gene_names, function(x) toString(x)))
allRes$gene_names  <- unlist(map(allRes$gene_names, function(x) sapply(strsplit(x, ","), paste, collapse=" ")))
#allResT <- allResT %>%
write.csv(allRes,str_glue("~/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/Cooltools_compartment_cluster/dcHiC_thres_{thres_decided}_clustered_k_means_cluster_{i}_GO_BP.csv"))}
##write.csv(allRes,str_glue("~/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/Cooltools_compartment_cluster/dcHiC_PC_clustered_cooltools_PCs_sample_wise_{nature}_k_means_cluster_{i}_GO_BP.csv"))}          
#write.csv(allRes,str_glue("~/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/dcHiC_Compartment_Cluster/clustered_dcHIC_PCs_{nature}_k_means_cluster_{i}_GO_BP.csv"))}
#write.csv(allRes,str_glue("~/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/Cooltools_compartment_cluster/clustered_cooltools_PCs_sample_wise_{nature}_k_means_cluster_{i}_GO_BP.csv"))}          
###

## load the dataset 
#data(geneList)
#library(package = affyLib, character.only = TRUE)

## the distribution of the adjusted p-values
#hist(geneList, 100)

## how many differentially expressed genes are:
#sum(topDiffGenes(geneList))

## build the topGOdata class 
#GOdata <- new("topGOdata",
              #ontology = "BP",
              #allGenes = geneList,
              #geneSel = topDiffGenes,
              #description = "GO analysis of ALL data: DE B-cell vs T-cell",
              #annot = annFUN.db,
              #affyLib = affyLib)

## display the GOdata object
#GOdata

##########################################################
## Examples on how to use the methods
##########################################################

## description of the experiment
#description(GOdata)

## obtain the genes that will be used in the analysis
#a <- genes(GOdata)
#str(a)
#numGenes(GOdata)

## obtain the score (p-value) of the genes
#selGenes <- names(geneList)[sample(1:length(geneList), 10)]
#gs <- geneScore(GOdata, whichGenes = selGenes)
#print(gs)

## if we want an unnamed vector containing all the feasible genes
#gs <- geneScore(GOdata, use.names = FALSE)
#str(gs)

## the list of significant genes
#sg <- sigGenes(GOdata)
#str(sg)
#numSigGenes(GOdata)

## to update the gene list 
#.geneList <- geneScore(GOdata, use.names = TRUE)
#GOdata ## more available genes
#GOdata <- updateGenes(GOdata, .geneList, topDiffGenes)
#GOdata ## the available genes are now the feasible genes

## the available GO terms (all the nodes in the graph)
#go <- usedGO(GOdata)
#length(go)

## to list the genes annotated to a set of specified GO terms
#sel.terms <- sample(go, 10)
#ann.genes <- genesInTerm(GOdata, sel.terms)
#str(ann.genes)

## the score for these genes
#ann.score <- scoresInTerm(GOdata, sel.terms)
#str(ann.score)

## to see the number of annotated genes
#num.ann.genes <- countGenesInTerm(GOdata)
#str(num.ann.genes)

## to summarise the statistics
#termStat(GOdata, sel.terms)
###
#listEnsembl()
#ensembl <- useEnsembl(biomart = "genes")
#datasets <- listDatasets(ensembl)
#head(datasets)
#listAttributes()

#regions <- data.frame(chr = query_df$chr, start = query_df$start, end = query_df$end)
#No need to use archived versions for now as ensembl is still using GRCz11 (May be changed to GRCz12Tu)
#ensembl <- useEnsembl(biomart = "genes", dataset = "drerio_gene_ensembl")
#filts<-listFilters(ensembl)
#attrs <- listAttributes(ensembl)
#results <- getBM(attributes = c("ensembl_gene_id","zfin_id_id"), 
      #filters = "chromosomal_region",
      #values = regions, 
      #mart = ensembl)

#P=Biological Process
#F=Molecular Function
#C=Cellular Component