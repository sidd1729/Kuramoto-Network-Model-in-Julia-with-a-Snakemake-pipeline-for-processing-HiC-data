#Gene_name to transcript choosing maximum length transcript
library(dplyr)
library(purrr)
Read_table <- read.table("~/Downloads/danRer10_refGene_gene.names.transcript.length.txt")
head(Read_table)
unique_genes <- unique(Read_table$V1)
unique_max_len_returner <- function(index){
  filtered_gene = unique_genes[index]
  filtered_table = Read_table %>% filter(V1 == filtered_gene)
  filtered_max_len = max(filtered_table$V2)
  return(filtered_max_len)
}

out_mapped <-map(seq(1,length(unique_genes),by = 1),unique_max_len_returner)
print(out_mapped)
gene_name_len_mapping <-data.frame(unique_genes = unique_genes, gene_lengths = unlist(out_mapped))
write.csv(gene_name_len_mapping,"~/Downloads/unique_danRer10_refGene_gene.names.transcript.length.txt",row.names = FALSE)

print(length(unique_genes) == length(out_mapped))