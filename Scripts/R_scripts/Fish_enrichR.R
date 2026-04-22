Sys.setenv(http_proxy = '172.16.2.252:3128')
Sys.setenv(https_proxy = '172.16.2.252:3128')

library(enrichR)
library(rtracklayer)
library(BRGenomics)
library(arrow)
library(ggplot2)
#library(purrr)
library(dplyr)
library(rbioapi)
library(stringr)
#data(genes790)

listEnrichrDbs()
listEnrichrSites()

#Setting enrichr site as Fishenrichr
setEnrichrSite("FishEnrichr")

#Databases to query in FishEnrichr
dbs <- c("GO_Molecular_Function_2018", "GO_Cellular_Component_2018",
                   "GO_Biological_Process_2018")

#Starting UCSC session to query Table Browser using R Tracklayer
mySession <- browserSession("UCSC")

#Setting genome as danRer11
genome(mySession) <- "danRer11"

#Loading in feather dataframe containing regions to intersect with track taken from table browser (not separate queries 
#as they were spaced by 10s to prevent excess requests if querying single regions and querying 1000 at a time just seemed to return the whole track )

cluster_no <- 7
feather_df <- read_feather("/home/ksslab/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/dcHiC_Compartment_Cluster/clustered_dcHIC_PCs_sample_wise_raw_k_means_df.feather")
feather_df <-feather_df %>% filter(cluster == cluster_no)
#Converting to GRange
snp_gr <-makeGRangesFromDataFrame(feather_df, seqnames.field = "chr", start.field = "start", end.field = "end")
for_query <- GenomicRanges::trim(snp_gr)

#Getting RefSeq curated track from UCSC Table Browser
tbl.rand <- getTable(ucscTableQuery(mySession,track = "NCBI RefSeq",table = "ncbiRefSeqCurated" ))

#Converting into GRange to intersect with query sequence
tbl.track_attempt <- makeGRangesFromDataFrame(tbl.rand,keep.extra.columns = TRUE,seqnames.field = "chrom", start.field = "txStart", end.field = "txEnd")
for_subject <- GenomicRanges::trim(tbl.track_attempt)

#Finding overlap of query with subject in subject and using subjectHits to index the hits/matching rows in the subject
overlap_obj<-findOverlaps(for_query,for_subject)
Cleared <- for_subject[subjectHits(overlap_obj)]

#Fishenrichr querying with extraction of gene name 
genes_to_query <- Cleared$name2
write(genes_to_query,"~/Downloads/gene_list_finally.txt")

#Submit the analysis request. from rbioapi tutorial
annots <- rba_panther_info(what = "datasets")
rba_panther_info(what = "organisms")
enriched <- rba_panther_enrich(
  genes = genes_to_query,
  organism = 7955,
  annot_dataset = "GO:0008150",
  cutoff = 0.05
)
result <- enriched$result
write.csv(result,str_glue("~/Downloads/diff_compartments__cluster_{cluster_no}_GO.csv"))
#> Performing PANTHER over-representation analysis (Fisher's exact test) on 15 genes from `organism 9606` against `GO:0008150` datasets.

# Note that we didn't supply the `test_type` parameter.
# In this case, the function will default to using Fisher's exact test # (i.e. `test_type = "FISHER"`).
# You may also use binomial test for the over-representation analysis # (i.e. `test_type = "BINOMIAL"`).

#EnrichR querying
df <- enrichr(
  genes_to_query,
  databases = dbs,
  background = NULL,
  include_overlap = FALSE,
  sleepTime = 1
)
write.csv(df$GO_Biological_Process_2018,str_glue("~/Downloads/diff_compartments__cluster_{cluster_no}_fishenrichR_GO_BP.csv"))
write.csv(df$GO_Cellular_Component_2018,str_glue("~/Downloads/diff_compartments__cluster_{cluster_no}_fishenrichR_GO_CC.csv"))
write.csv(df$GO_Molecular_Function_2018,str_glue("~/Downloads/diff_compartments__cluster_{cluster_no}_fishenrichR_GO_MF.csv"))

#Plotting results of any one of the databases queried
png(str_glue("~/Downloads/diff_compartments_cluster_{cluster_no}_fishenrichR_GO_BP.png"))
plotEnrich(
  df$GO_Biological_Process_2018,
  showTerms = 20,
  numChar = 60,
  y = "Count",
  orderBy = "P.value",
  xlab = NULL,
  ylab = NULL,
  title = str_glue("Cluster {cluster_no} BP enrichment ")
)
dev.off()
png(str_glue("~/Downloads/diff_compartments_cluster_{cluster_no}_fishenrichR_GO_CC.png"))
plotEnrich(
  df$GO_Cellular_Component_2018,
  showTerms = 20,
  numChar = 60,
  y = "Count",
  orderBy = "P.value",
  xlab = NULL,
  ylab = NULL,
  title = str_glue("Cluster {cluster_no} CC enrichment ")
)
dev.off()
png(str_glue("~/Downloads/diff_compartments_cluster_{cluster_no}_fishenrichR_GO_MF.png"))
plotEnrich(
  df$GO_Molecular_Function_2018,
  showTerms = 20,
  numChar = 60,
  y = "Count",
  orderBy = "P.value",
  xlab = NULL,
  ylab = NULL,
  title = str_glue("Cluster {cluster_no} MF enrichment")
)
dev.off()


#snp_gr <- renameSeqlevels(snp_gr, paste0("chr", seqlevels(snp_gr)))
#for_query <- GRangesForUCSCGenome("danRer11", chrom = snp_gr@seqnames, ranges = snp_gr@ranges)
#for_query <- GenomicRanges::trim(for_query)
#Converting Dataframe to GRanges object to query the Tablebrowser
#rand.tss2.grange <- GRanges("chr5",IRanges(40000000,41000000))
#rand.tss.grange <- GRanges("chr19",IRanges(15000000,30000000))
#gr_list <- GRangesList("gr1" = rand.tss.grange, "gr2" = rand.tss2.grange)
#range_len <- dim(feather_df)[1]
#seq_map <- seq(1,range_len,by=1)
#chr_vec <- feather_df$chr
#start_vec <- feather_df$start
#end_vec <- feather_df$end
#GRange_list_hopefully <- map(seq_map,Gene_name_extractor)
#Gene_name_extractor <- function(x){
  #grange_ex <- GRanges(chr_vec[x],IRanges(start_vec[x],end_vec[x]))
  #tbl.rand <- getTable(ucscTableQuery(mySession,track = "NCBI RefSeq",range = grange_ex, table = "ncbiRefSeqCurated" ))
  #Sys.sleep(10)
  #return (tbl.rand$name2)
#}
#snp_gr <-makeGRangesFromDataFrame(feather_df, seqnames.field = "chr", start.field = "start", end.field = "end")
#snp_gr <- renameSeqlevels(snp_gr, paste0("chr", seqlevels(snp_gr)))
#for_query <- GRangesForUCSCGenome("danRer11", chrom = snp_gr@seqnames, ranges = snp_gr@ranges)
#for_query <- GenomicRanges::trim(for_query)

#randmerge.tss.grange <-pc(rand.tss.grange,rand.tss2.grange)

#if (range_len > 1000){
#tbl.rand <- getTable(ucscTableQuery(mySession,track = "NCBI RefSeq",range = for_query[1:1000] ,table = "ncbiRefSeqCurated" ))
#tbl.rand2 <- getTable(ucscTableQuery(mySession,track = "NCBI RefSeq",range = for_query[1001:range_len] ,table = "ncbiRefSeqCurated" ))}else{
  #tbl.rand <- getTable(ucscTableQuery(mySession,track = "NCBI RefSeq",table = "ncbiRefSeqCurated" ))
  #tbl.track <- rtracklayer::track(ucscTableQuery(mySession,track = "NCBI RefSeq",table = "ncbiRefSeqCurated" ))#,asRangedData=FALSE)
#}
#tbl.track_attempt <- makeGRangesFromDataFrame(tbl.rand,keep.extra.columns = TRUE,seqnames.field = "chrom", start.field = "txStart", end.field = "txEnd")
#for_subject <- GRangesForUCSCGenome("danRer11", chrom = tbl.track_attempt@seqnames, ranges = tbl.track_attempt@ranges)
#for_subject <- trim(for_subject)
#Tested <-disjoin(c(for_query,tbl.track_attempt), with.revmap = TRUE)
#revmap <- Tested$revmap
#r_scores <- extractList(mcols(c(for_query,tbl.track_attempt))$name2,revmap)
#length(for_query) == range_len
#overlap_obj<-findOverlaps(for_query,tbl.track_attempt)
#Cleared <- tbl.track_attempt[subjectHits(overlap_obj)]
#for (i in seq(0,range_len%/%1000,by=1)){
  #print(i)
  #print(1+(1000*i))
  #print(1000*(i+1))
  #if (range_len >= 1000*(i+1)){
    #print("Yeah")
    #print(i)
    #print(1+(1000*i))
    #print(1000*(i+1))
  #tbl.rand <- getTable(ucscTableQuery(mySession,track = "NCBI RefSeq",range = for_query[(1+(1000*i)):(1000*(i+1))] ,table = "ncbiRefSeqCurated" ))
  #print("Yeah 1000 done")
  #print(i)
  #print(1+(1000*i))
  #print(1000*(i+1))}
  #else{
    #print("Nah")
    #print(i)
    #print(1+(1000*i))
    #print(1000*(i+1))
    #print((1+(1000*i)):range_len)
    #tbl.rand <- getTable(ucscTableQuery(mySession,track = "NCBI RefSeq",range = for_query[(1+(1000*i)):range_len] ,table = "ncbiRefSeqCurated" ))
    #print(i)
    #print(1+(1000*i))
    #print(1000*(i+1))
  #print("Nah, less than 1000")}
#}
library(maotai)

kmeanspp()
#tbl.track <- GRanges(ucscTableQuery(mySession,track = "NCBI RefSeq", table = "ncbiRefSeqCurated" ))
genes_2 <- tbl.rand$name2
genes <- c("dtnbp1b", "xkr8.3", "col28a1a", "tmem268", "tshz1", "tbc1d5", "vars2", "csnk2b",
           "si:ch211-30b16.2", "osbpl3a", "nelfe", "hoxa3a", "npy", "ptpn23b", "hdac1",
           "lpcat1", "prrc2a", "npr1a", "adcy2b", "npvf", "jazf1a", "ctps1a", "prf1.2",
           "zgc:64022", "mrpl3", "mylipa", "LOC108180149", "ndufs6", "med10", "snx10a",
           "nudcd1", "mir219-3", "LOC101882396", "znf516", "gdf6b", "glcci1", "hoxa4a",
           "mir196a-2", "top2b", "ngly1", "fndc5b", "ddah2", "hibadha", "grb10a", "znrd1",
           "atp9b", "tac1", "rab5aa", "srd5a1", "adnp2b", "txnl4a", "txnipa", "hoxa1a",
           "pard6gb", "si:dkey-119m17.2", "hoxa9a", "cndp2", "satb1a", "stmn1a",
           "si:ch211-206a7.2", "si:dkey-261i16.5", "evx1", "slc39a7", "si:dkey-81h8.1",
           "ppt2", "bloc1s4", "ptp4a3", "polr3gla", "si:dkey-202e17.1", "fbxl6", "chtopa",
           "ldlrap1a", "smarcc1b", "gmeb1", "irgf1", "cratb", "nkiras1", "thap7", "rps18",
           "irx4b", "pkhd1l1", "tmem57a", "hoxa11a", "skiv2l", "lin28a", "nutf2l", "ints3",
           "si:dkey-208k4.2", "fam221a", "dazl", "txlna", "rims3", "mgat1b",
           "si:ch211-244a23.1", "dync1li1", "edn2", "si:ch211-195b13.1", "serinc2",
           "zgc:100906", "gtf2h4", "ccdc126", "nsun2", "eif3i", "zgc:122991", "zbtb12.1",
           "sall3b", "zbtb12.2", "asns", "pex11b", "igf2bp3", "ddx39b", "cited4a", "s100v2",
           "pdik1l", "rnu11", "ppp1r8a", "rab42a", "mir7145", "cbx3a", "clic1", "srfbp1",
           "tax1bp1a", "syt11a", "ppp1r11", "nfe2l3", "marcksl1b", "rpl15", "mpp6b",
           "sf3a3", "slc52a2", "sacm1la", "hoxa5a", "jarid2b", "mtx1b", "hoxa13a",
           "mrpl36", "sdhaf3", "p3h4", "sox4a", "trhrb", "tnfa", "irx1b", "rit1", "tra2a",
           "atat1", "eppk1", "eya3", "thbs3b", "cspg5b", "si:dkey-222b8.1", "oxnad1",
           "thrb", "nr1d2b", "snapin", "e2f3", "nfatc1", "rpa3", "si:dkeyp-92c9.2",
           "kat2b", "eif3eb", "ube2e1")

df <- enrichr(
  unlisted_trial,
  databases = dbs,
  background = NULL,
  include_overlap = FALSE,
  sleepTime = 1
)
plotEnrich(
  df$GO_Molecular_Function_2018,
  showTerms = 20,
  numChar = 40,
  y = "Count",
  orderBy = "P.value",
  xlab = NULL,
  ylab = NULL,
  title = NULL
)
test_attempt <- unique(unlist(r_scores))
saveRDS(GRange_list_hopefully,"~/Downloads/gene_list_finally.rds")
trial<-readRDS("~/Downloads/gene_list_finally.rds")
unlisted_trial <- unlist(trial)
test <- read_feather("~/Siddharth/Program_Directory/Kuramoto_Processing_Pipeline/results/dcHiC_Compartment_Cluster/clustered_dcHIC_PCs_sample_wise_raw_k_means_centers_df.feather")