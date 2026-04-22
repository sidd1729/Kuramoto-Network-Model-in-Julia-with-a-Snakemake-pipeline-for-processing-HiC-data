#First converting the gene IDs to gene names
library(tidyr)
library(stringr)
library(org.Hs.eg.db)
library(annotables)
library(biomaRt)
library(rentrez)
library(furrr)
library(purrr)
library(arrow)
library(data.table)
library(dplyr)
library(dbplyr)
#Reading in our TSV file of counts
GM12878_counts <- tibble(read.delim("~/Siddharth/Program_Directory/Data/RNA_Seq_Data/gene_quant/ENCFF513SFV.tsv"))
GM12878_counts <- GM12878_counts %>% dplyr::filter(TPM != 0)
GM12878_counts$gene_id <- gsub('\\..*','',GM12878_counts$gene_id)
#GM12878_counts_ID <- gsub('\\..*','',GM12878_counts_ID)
GM12878_counts %>% separate_longer_delim(gene_id,delim ="ENS")
GM12878_counts %>% separate_longer_delim(gene_id,delim ="\\d+")
ENS_IDs <- GM12878_counts %>% dplyr::filter(str_detect(gene_id, "ENS.+"))
ENT_IDs <- GM12878_counts %>% dplyr::filter(str_detect(gene_id, "ENS.+",negate = T))
Truly_ENT_IDs <- ENT_IDs %>% dplyr::filter(str_detect(gene_id,"[a-zA-Z]+",negate = T))
Random_other_IDs <- ENT_IDs %>% dplyr::filter(str_detect(gene_id,"[a-zA-Z]+"))
ENS_mapped <- data.frame(mapIds(org.Hs.eg.db,keys = ENS_IDs$gene_id, keytype = "ENSEMBL", column = "SYMBOL"))
ENT_mapped <- data.frame(mapIds(org.Hs.eg.db,keys = ENT_IDs$gene_id, keytype = "ENTREZID",column = "SYMBOL"))
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
gene_pos <- getBM(attributes=c('ensembl_gene_id','external_gene_name','chromosome_name','start_position','end_position'),filters = 'ensembl_gene_id',values = ENS_IDs$gene_id,mart =ensembl)
gene_pos_entr <- getBM(attributes = c('entrezgene_id','external_gene_name','chromosome_name','start_position','end_position'),filters ='entrezgene_id', values = ENT_IDs$gene_id, mart = ensembl)

#Unmapped Ensembl IDs (Cartesian/Set Difference between total and mapped sets:
umapped_ids  <- setdiff(ENS_IDs$gene_id, gene_pos$ensembl_gene_id)
write.table(umapped_ids,"~/Siddharth/Program_Directory/Data/RNA_Seq_Data/gene_quant/unmapped_ensembl.txt",row.names = FALSE,col.names = FALSE,quote = FALSE)
entrez=c("673","837")
goids = getBM(attributes=c('entrezgene_id','go_id'), filters='entrezgene_id', values=entrez, mart=ensembl)
head(goids)
filter(df, id == 1)
colnames(GM12878_counts)
mapIds()

write.table(ENT_IDs$gene_id,"~/Siddharth/Program_Directory/Data/RNA_Seq_Data/gene_quant/entrez.txt",row.names = FALSE,col.names = FALSE,quote = FALSE)
write.table(ENS_IDs$gene_id,"~/Siddharth/Program_Directory/Data/RNA_Seq_Data/gene_quant/ensembl.txt",row.names = FALSE,col.names = FALSE,quote = FALSE)
# if(interactive()){
#   mart <- useMart(biomart = "ensembl", 
#                      dataset = "hsapiens_gene_ensembl")
#   
#   getBM(attributes = c("affy_hg_u95av2", "hgnc_symbol", "chromosome_name", "band"),
#         filters    = "affy_hg_u95av2",
#         values     = c("1939_at","1503_at","1454_at"), 
#         mart       = mart)
# }
 gene_link <- entrez_link(dbfrom = "gene",id = 306, db ="nuccore")
 esum <-entrez_summary(db = "gene",id = 306)
 
 Entrez_pos_retriever <- function(gene_id_val){
 gene_link <- entrez_link(dbfrom = "gene",id = gene_id_val, db ="nuccore")
 esum <-entrez_summary(db = "gene",id = gene_id_val)
 return (list(gene_id_val,esum$name,esum$genomicinfo$chrloc, esum$genomicinfo$chrstart,esum$genomicinfo$chrstop))}

plan(multisession, workers = 30,gc=TRUE)
Coord_mapped <-map(Truly_ENT_IDs$gene_id,Entrez_pos_retriever)
colnames(gene_pos)[1] <- "gene_id"

first_test <- right_join(ENS_IDs,gene_pos,by = join_by(gene_id))

Coord_df <- data.frame(Coord_mapped)

tied_together <-data.table::rbindlist(list(gene_pos, Coord_df),use.names=FALSE)
overall_params_df <- GM12878_counts[c("gene_id","FPKM","TPM")]
colnames(overall_params_df)[1] <- "ensembl_gene_id"
express_pos_df <- right_join(overall_params_df,tied_together,by = join_by(ensembl_gene_id))
write_feather(express_pos_df,"~/Siddharth/Program_Directory/Data/RNA_Seq_Data/gene_quant/plottable_params_GM12878.feather")
