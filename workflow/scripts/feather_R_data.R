library(arrow)
library(dplyr)
library(stringr)
options(scipen=999)
sample_wise_df_2 <- read_feather("results/k_50_table/k_50_sample_wise_df.feather")
df_to_hmap_2 <- read_feather("results/k_50_table/filtered_non_clustered_k_50_sample_wise_df.feather")
#df_to_hmap_2$chromosome <- NA 
#df_to_hmap_2$region <- NA
df_to_hmap_2 <- df_to_hmap_2 %>% mutate(Chromosome = case_when( as.numeric(bin_no) <= 1188 ~ 1, 
as.numeric(bin_no) > 1188 & as.numeric(bin_no) <= 2377 ~ 2, as.numeric(bin_no) > 2377 & as.numeric(bin_no) <= 3626 ~ 3, 
as.numeric(bin_no) > 3626 & as.numeric(bin_no) <= 5184 ~ 4, as.numeric(bin_no) > 5184 & as.numeric(bin_no) <= 6631 ~ 5,
as.numeric(bin_no) > 6631 & as.numeric(bin_no) <= 7833 ~ 6, as.numeric(bin_no) > 7833 & as.numeric(bin_no) <= 9315 ~ 7,
as.numeric(bin_no) > 9315 & as.numeric(bin_no) <= 10398 ~ 8, as.numeric(bin_no) > 10398 & as.numeric(bin_no) <= 11524 ~ 9,
as.numeric(bin_no) > 11524 & as.numeric(bin_no) <= 12429 ~ 10, as.numeric(bin_no) > 12429 & as.numeric(bin_no) <= 13335 ~ 11,
as.numeric(bin_no) > 13335 & as.numeric(bin_no) <= 14315 ~ 12, as.numeric(bin_no) > 14315 & as.numeric(bin_no) <= 15355 ~ 13,
as.numeric(bin_no) > 15355 & as.numeric(bin_no) <= 16405 ~ 14, as.numeric(bin_no) > 16405 & as.numeric(bin_no) <= 17362 ~ 15,
as.numeric(bin_no) > 17362 & as.numeric(bin_no) <= 18464 ~ 16, as.numeric(bin_no) > 18464 & as.numeric(bin_no) <= 19530 ~ 17,
as.numeric(bin_no) > 19530 & as.numeric(bin_no) <= 20547 ~ 18, as.numeric(bin_no) > 20547 & as.numeric(bin_no) <= 21512 ~ 19,
as.numeric(bin_no) > 21512 & as.numeric(bin_no) <= 22613 ~ 20, as.numeric(bin_no) > 22613 & as.numeric(bin_no) <= 23528 ~ 21,
as.numeric(bin_no) > 23528 & as.numeric(bin_no) <= 24307 ~ 22, as.numeric(bin_no) > 24307 & as.numeric(bin_no) <= 25228 ~ 23,
as.numeric(bin_no) > 25228 & as.numeric(bin_no) <= 26068 ~ 24, as.numeric(bin_no) > 26068 & as.numeric(bin_no) <= 26815 ~ 25))

df_to_hmap_2 <- df_to_hmap_2 %>% mutate(Region = case_when( as.numeric(bin_no) <= 1188 ~ str_glue("{(as.numeric(bin_no)-1)*50000}-{(as.numeric(bin_no))*50000}"), 
as.numeric(bin_no) > 1188 & as.numeric(bin_no) <= 2377 ~ str_glue("{(as.numeric(bin_no)-1-1188)*50000}-{(as.numeric(bin_no)-1188)*50000}"), as.numeric(bin_no) > 2377 & as.numeric(bin_no) <= 3626 ~ str_glue("{(as.numeric(bin_no)-1-2377)*50000}-{(as.numeric(bin_no)-2377)*50000}"), 
as.numeric(bin_no) > 3626 & as.numeric(bin_no) <= 5184 ~ str_glue("{(as.numeric(bin_no)-1-3626)*50000}-{(as.numeric(bin_no)-3626)*50000}"), as.numeric(bin_no) > 5184 & as.numeric(bin_no) <= 6631 ~ str_glue("{(as.numeric(bin_no)-1-5184)*50000}-{(as.numeric(bin_no)-5184)*50000}"),
as.numeric(bin_no) > 6631 & as.numeric(bin_no) <= 7833 ~ str_glue("{(as.numeric(bin_no)-1-6631)*50000}-{(as.numeric(bin_no)-6631)*50000}"), as.numeric(bin_no) > 7833 & as.numeric(bin_no) <= 9315 ~ str_glue("{(as.numeric(bin_no)-1-7833)*50000}-{(as.numeric(bin_no)-7833)*50000}"),
as.numeric(bin_no) > 9315 & as.numeric(bin_no) <= 10398 ~ str_glue("{(as.numeric(bin_no)-1-9315)*50000}-{(as.numeric(bin_no)-9315)*50000}"), as.numeric(bin_no) > 10398 & as.numeric(bin_no) <= 11524 ~ str_glue("{(as.numeric(bin_no)-1-10398)*50000}-{(as.numeric(bin_no)-10398)*50000}"),
as.numeric(bin_no) > 11524 & as.numeric(bin_no) <= 12429 ~ str_glue("{(as.numeric(bin_no)-1-11524)*50000}-{(as.numeric(bin_no)-11524)*50000}"), as.numeric(bin_no) > 12429 & as.numeric(bin_no) <= 13335 ~ str_glue("{(as.numeric(bin_no)-1-12429)*50000}-{(as.numeric(bin_no)-12429)*50000}"),
as.numeric(bin_no) > 13335 & as.numeric(bin_no) <= 14315 ~ str_glue("{(as.numeric(bin_no)-1-13335)*50000}-{(as.numeric(bin_no)-13335)*50000}"), as.numeric(bin_no) > 14315 & as.numeric(bin_no) <= 15355 ~ str_glue("{(as.numeric(bin_no)-1-14315)*50000}-{(as.numeric(bin_no)-14315)*50000}"),
as.numeric(bin_no) > 15355 & as.numeric(bin_no) <= 16405 ~ str_glue("{(as.numeric(bin_no)-1-15355)*50000}-{(as.numeric(bin_no)-15355)*50000}"), as.numeric(bin_no) > 16405 & as.numeric(bin_no) <= 17362 ~ str_glue("{(as.numeric(bin_no)-1-16405)*50000}-{(as.numeric(bin_no)-16405)*50000}"),
as.numeric(bin_no) > 17362 & as.numeric(bin_no) <= 18464 ~ str_glue("{(as.numeric(bin_no)-1-17362)*50000}-{(as.numeric(bin_no)-17362)*50000}"), as.numeric(bin_no) > 18464 & as.numeric(bin_no) <= 19530 ~ str_glue("{(as.numeric(bin_no)-1-18464)*50000}-{(as.numeric(bin_no)-18464)*50000}"),
as.numeric(bin_no) > 19530 & as.numeric(bin_no) <= 20547 ~ str_glue("{(as.numeric(bin_no)-1-19530)*50000}-{(as.numeric(bin_no)-19530)*50000}"), as.numeric(bin_no) > 20547 & as.numeric(bin_no) <= 21512 ~ str_glue("{(as.numeric(bin_no)-1- 20547)*50000}-{(as.numeric(bin_no)-20547)*50000}"),
as.numeric(bin_no) > 21512 & as.numeric(bin_no) <= 22613 ~ str_glue("{(as.numeric(bin_no)-1-21512)*50000}-{(as.numeric(bin_no)-21512)*50000}"), as.numeric(bin_no) > 22613 & as.numeric(bin_no) <= 23528 ~ str_glue("{(as.numeric(bin_no)-1-22613)*50000}-{(as.numeric(bin_no)-22613)*50000}"),
as.numeric(bin_no) > 23528 & as.numeric(bin_no) <= 24307 ~ str_glue("{(as.numeric(bin_no)-1-23528)*50000}-{(as.numeric(bin_no) - 23528)*50000}"), as.numeric(bin_no) > 24307 & as.numeric(bin_no) <= 25228 ~ str_glue("{(as.numeric(bin_no)-1 - 24307)*50000}-{(as.numeric(bin_no)-24307)*50000}"),
as.numeric(bin_no) > 25228 & as.numeric(bin_no) <= 26068 ~ str_glue("{(as.numeric(bin_no)-1-25228)*50000}-{(as.numeric(bin_no)-25228)*50000}"), as.numeric(bin_no) > 26068 & as.numeric(bin_no) <= 26815 ~ str_glue("{(as.numeric(bin_no)-1-26068)*50000}-{(as.numeric(bin_no)-26068)*50000}")))

df_to_hmap_2 <- df_to_hmap_2 %>% mutate(Chromosome_Region = case_when( as.numeric(bin_no) <= 1188 ~ str_glue("chr{Chromosome}:{Region}"), 
as.numeric(bin_no) > 1188 & as.numeric(bin_no) <= 2377 ~ str_glue("chr{Chromosome}:{Region}"), as.numeric(bin_no) > 2377 & as.numeric(bin_no) <= 3626 ~ str_glue("chr{Chromosome}:{Region}"), 
as.numeric(bin_no) > 3626 & as.numeric(bin_no) <= 5184 ~ str_glue("chr{Chromosome}:{Region}"), as.numeric(bin_no) > 5184 & as.numeric(bin_no) <= 6631 ~ str_glue("chr{Chromosome}:{Region}"),
as.numeric(bin_no) > 6631 & as.numeric(bin_no) <= 7833 ~ str_glue("chr{Chromosome}:{Region}"), as.numeric(bin_no) > 7833 & as.numeric(bin_no) <= 9315 ~ str_glue("chr{Chromosome}:{Region}"),
as.numeric(bin_no) > 9315 & as.numeric(bin_no) <= 10398 ~ str_glue("chr{Chromosome}:{Region}"), as.numeric(bin_no) > 10398 & as.numeric(bin_no) <= 11524 ~ str_glue("chr{Chromosome}:{Region}"),
as.numeric(bin_no) > 11524 & as.numeric(bin_no) <= 12429 ~ str_glue("chr{Chromosome}:{Region}"), as.numeric(bin_no) > 12429 & as.numeric(bin_no) <= 13335 ~ str_glue("chr{Chromosome}:{Region}"),
as.numeric(bin_no) > 13335 & as.numeric(bin_no) <= 14315 ~ str_glue("chr{Chromosome}:{Region}"), as.numeric(bin_no) > 14315 & as.numeric(bin_no) <= 15355 ~ str_glue("chr{Chromosome}:{Region}"),
as.numeric(bin_no) > 15355 & as.numeric(bin_no) <= 16405 ~ str_glue("chr{Chromosome}:{Region}"), as.numeric(bin_no) > 16405 & as.numeric(bin_no) <= 17362 ~ str_glue("chr{Chromosome}:{Region}"),
as.numeric(bin_no) > 17362 & as.numeric(bin_no) <= 18464 ~ str_glue("chr{Chromosome}:{Region}"), as.numeric(bin_no) > 18464 & as.numeric(bin_no) <= 19530 ~ str_glue("chr{Chromosome}:{Region}"),
as.numeric(bin_no) > 19530 & as.numeric(bin_no) <= 20547 ~ str_glue("chr{Chromosome}:{Region}"), as.numeric(bin_no) > 20547 & as.numeric(bin_no) <= 21512 ~ str_glue("chr{Chromosome}:{Region}"),
as.numeric(bin_no) > 21512 & as.numeric(bin_no) <= 22613 ~ str_glue("chr{Chromosome}:{Region}"), as.numeric(bin_no) > 22613 & as.numeric(bin_no) <= 23528 ~ str_glue("chr{Chromosome}:{Region}"),
as.numeric(bin_no) > 23528 & as.numeric(bin_no) <= 24307 ~ str_glue("chr{Chromosome}:{Region}"), as.numeric(bin_no) > 24307 & as.numeric(bin_no) <= 25228 ~ str_glue("chr{Chromosome}:{Region}"),
as.numeric(bin_no) > 25228 & as.numeric(bin_no) <= 26068 ~ str_glue("chr{Chromosome}:{Region}"), as.numeric(bin_no) > 26068 & as.numeric(bin_no) <= 26815 ~ str_glue("chr{Chromosome}:{Region}")))

#df_to_hmap_2$Chromosome <- NULL 
#df_to_hmap_2$Region <- NULL 

sample_wise_df_2$bin_no <- rownames(sample_wise_df_2)
sample_wise_df_2 <- sample_wise_df_2 %>% mutate(Chromosome = case_when( as.numeric(bin_no) <= 1188 ~ 1, 
as.numeric(bin_no) > 1188 & as.numeric(bin_no) <= 2377 ~ 2, as.numeric(bin_no) > 2377 & as.numeric(bin_no) <= 3626 ~ 3, 
as.numeric(bin_no) > 3626 & as.numeric(bin_no) <= 5184 ~ 4, as.numeric(bin_no) > 5184 & as.numeric(bin_no) <= 6631 ~ 5,
as.numeric(bin_no) > 6631 & as.numeric(bin_no) <= 7833 ~ 6, as.numeric(bin_no) > 7833 & as.numeric(bin_no) <= 9315 ~ 7,
as.numeric(bin_no) > 9315 & as.numeric(bin_no) <= 10398 ~ 8, as.numeric(bin_no) > 10398 & as.numeric(bin_no) <= 11524 ~ 9,
as.numeric(bin_no) > 11524 & as.numeric(bin_no) <= 12429 ~ 10, as.numeric(bin_no) > 12429 & as.numeric(bin_no) <= 13335 ~ 11,
as.numeric(bin_no) > 13335 & as.numeric(bin_no) <= 14315 ~ 12, as.numeric(bin_no) > 14315 & as.numeric(bin_no) <= 15355 ~ 13,
as.numeric(bin_no) > 15355 & as.numeric(bin_no) <= 16405 ~ 14, as.numeric(bin_no) > 16405 & as.numeric(bin_no) <= 17362 ~ 15,
as.numeric(bin_no) > 17362 & as.numeric(bin_no) <= 18464 ~ 16, as.numeric(bin_no) > 18464 & as.numeric(bin_no) <= 19530 ~ 17,
as.numeric(bin_no) > 19530 & as.numeric(bin_no) <= 20547 ~ 18, as.numeric(bin_no) > 20547 & as.numeric(bin_no) <= 21512 ~ 19,
as.numeric(bin_no) > 21512 & as.numeric(bin_no) <= 22613 ~ 20, as.numeric(bin_no) > 22613 & as.numeric(bin_no) <= 23528 ~ 21,
as.numeric(bin_no) > 23528 & as.numeric(bin_no) <= 24307 ~ 22, as.numeric(bin_no) > 24307 & as.numeric(bin_no) <= 25228 ~ 23,
as.numeric(bin_no) > 25228 & as.numeric(bin_no) <= 26068 ~ 24, as.numeric(bin_no) > 26068 & as.numeric(bin_no) <= 26815 ~ 25))

sample_wise_df_2 <- sample_wise_df_2 %>% mutate(Region = case_when( as.numeric(bin_no) <= 1188 ~ str_glue("{(as.numeric(bin_no)-1)*50000}-{(as.numeric(bin_no))*50000}"), 
as.numeric(bin_no) > 1188 & as.numeric(bin_no) <= 2377 ~ str_glue("{(as.numeric(bin_no)-1-1188)*50000}-{(as.numeric(bin_no)-1188)*50000}"), as.numeric(bin_no) > 2377 & as.numeric(bin_no) <= 3626 ~ str_glue("{(as.numeric(bin_no)-1-2377)*50000}-{(as.numeric(bin_no)-2377)*50000}"), 
as.numeric(bin_no) > 3626 & as.numeric(bin_no) <= 5184 ~ str_glue("{(as.numeric(bin_no)-1-3626)*50000}-{(as.numeric(bin_no)-3626)*50000}"), as.numeric(bin_no) > 5184 & as.numeric(bin_no) <= 6631 ~ str_glue("{(as.numeric(bin_no)-1-5184)*50000}-{(as.numeric(bin_no)-5184)*50000}"),
as.numeric(bin_no) > 6631 & as.numeric(bin_no) <= 7833 ~ str_glue("{(as.numeric(bin_no)-1-6631)*50000}-{(as.numeric(bin_no)-6631)*50000}"), as.numeric(bin_no) > 7833 & as.numeric(bin_no) <= 9315 ~ str_glue("{(as.numeric(bin_no)-1-7833)*50000}-{(as.numeric(bin_no)-7833)*50000}"),
as.numeric(bin_no) > 9315 & as.numeric(bin_no) <= 10398 ~ str_glue("{(as.numeric(bin_no)-1-9315)*50000}-{(as.numeric(bin_no)-9315)*50000}"), as.numeric(bin_no) > 10398 & as.numeric(bin_no) <= 11524 ~ str_glue("{(as.numeric(bin_no)-1-10398)*50000}-{(as.numeric(bin_no)-10398)*50000}"),
as.numeric(bin_no) > 11524 & as.numeric(bin_no) <= 12429 ~ str_glue("{(as.numeric(bin_no)-1-11524)*50000}-{(as.numeric(bin_no)-11524)*50000}"), as.numeric(bin_no) > 12429 & as.numeric(bin_no) <= 13335 ~ str_glue("{(as.numeric(bin_no)-1-12429)*50000}-{(as.numeric(bin_no)-12429)*50000}"),
as.numeric(bin_no) > 13335 & as.numeric(bin_no) <= 14315 ~ str_glue("{(as.numeric(bin_no)-1-13335)*50000}-{(as.numeric(bin_no)-13335)*50000}"), as.numeric(bin_no) > 14315 & as.numeric(bin_no) <= 15355 ~ str_glue("{(as.numeric(bin_no)-1-14315)*50000}-{(as.numeric(bin_no)-14315)*50000}"),
as.numeric(bin_no) > 15355 & as.numeric(bin_no) <= 16405 ~ str_glue("{(as.numeric(bin_no)-1-15355)*50000}-{(as.numeric(bin_no)-15355)*50000}"), as.numeric(bin_no) > 16405 & as.numeric(bin_no) <= 17362 ~ str_glue("{(as.numeric(bin_no)-1-16405)*50000}-{(as.numeric(bin_no)-16405)*50000}"),
as.numeric(bin_no) > 17362 & as.numeric(bin_no) <= 18464 ~ str_glue("{(as.numeric(bin_no)-1-17362)*50000}-{(as.numeric(bin_no)-17362)*50000}"), as.numeric(bin_no) > 18464 & as.numeric(bin_no) <= 19530 ~ str_glue("{(as.numeric(bin_no)-1-18464)*50000}-{(as.numeric(bin_no)-18464)*50000}"),
as.numeric(bin_no) > 19530 & as.numeric(bin_no) <= 20547 ~ str_glue("{(as.numeric(bin_no)-1-19530)*50000}-{(as.numeric(bin_no)-19530)*50000}"), as.numeric(bin_no) > 20547 & as.numeric(bin_no) <= 21512 ~ str_glue("{(as.numeric(bin_no)-1- 20547)*50000}-{(as.numeric(bin_no)-20547)*50000}"),
as.numeric(bin_no) > 21512 & as.numeric(bin_no) <= 22613 ~ str_glue("{(as.numeric(bin_no)-1-21512)*50000}-{(as.numeric(bin_no)-21512)*50000}"), as.numeric(bin_no) > 22613 & as.numeric(bin_no) <= 23528 ~ str_glue("{(as.numeric(bin_no)-1-22613)*50000}-{(as.numeric(bin_no)-22613)*50000}"),
as.numeric(bin_no) > 23528 & as.numeric(bin_no) <= 24307 ~ str_glue("{(as.numeric(bin_no)-1-23528)*50000}-{(as.numeric(bin_no) - 23528)*50000}"), as.numeric(bin_no) > 24307 & as.numeric(bin_no) <= 25228 ~ str_glue("{(as.numeric(bin_no)-1 - 24307)*50000}-{(as.numeric(bin_no)-24307)*50000}"),
as.numeric(bin_no) > 25228 & as.numeric(bin_no) <= 26068 ~ str_glue("{(as.numeric(bin_no)-1-25228)*50000}-{(as.numeric(bin_no)-25228)*50000}"), as.numeric(bin_no) > 26068 & as.numeric(bin_no) <= 26815 ~ str_glue("{(as.numeric(bin_no)-1-26068)*50000}-{(as.numeric(bin_no)-26068)*50000}")))

sample_wise_df_2 <- sample_wise_df_2 %>% mutate(Chromosome_Region = case_when( as.numeric(bin_no) <= 1188 ~ str_glue("chr{Chromosome}:{Region}"), 
as.numeric(bin_no) > 1188 & as.numeric(bin_no) <= 2377 ~ str_glue("chr{Chromosome}:{Region}"), as.numeric(bin_no) > 2377 & as.numeric(bin_no) <= 3626 ~ str_glue("chr{Chromosome}:{Region}"), 
as.numeric(bin_no) > 3626 & as.numeric(bin_no) <= 5184 ~ str_glue("chr{Chromosome}:{Region}"), as.numeric(bin_no) > 5184 & as.numeric(bin_no) <= 6631 ~ str_glue("chr{Chromosome}:{Region}"),
as.numeric(bin_no) > 6631 & as.numeric(bin_no) <= 7833 ~ str_glue("chr{Chromosome}:{Region}"), as.numeric(bin_no) > 7833 & as.numeric(bin_no) <= 9315 ~ str_glue("chr{Chromosome}:{Region}"),
as.numeric(bin_no) > 9315 & as.numeric(bin_no) <= 10398 ~ str_glue("chr{Chromosome}:{Region}"), as.numeric(bin_no) > 10398 & as.numeric(bin_no) <= 11524 ~ str_glue("chr{Chromosome}:{Region}"),
as.numeric(bin_no) > 11524 & as.numeric(bin_no) <= 12429 ~ str_glue("chr{Chromosome}:{Region}"), as.numeric(bin_no) > 12429 & as.numeric(bin_no) <= 13335 ~ str_glue("chr{Chromosome}:{Region}"),
as.numeric(bin_no) > 13335 & as.numeric(bin_no) <= 14315 ~ str_glue("chr{Chromosome}:{Region}"), as.numeric(bin_no) > 14315 & as.numeric(bin_no) <= 15355 ~ str_glue("chr{Chromosome}:{Region}"),
as.numeric(bin_no) > 15355 & as.numeric(bin_no) <= 16405 ~ str_glue("chr{Chromosome}:{Region}"), as.numeric(bin_no) > 16405 & as.numeric(bin_no) <= 17362 ~ str_glue("chr{Chromosome}:{Region}"),
as.numeric(bin_no) > 17362 & as.numeric(bin_no) <= 18464 ~ str_glue("chr{Chromosome}:{Region}"), as.numeric(bin_no) > 18464 & as.numeric(bin_no) <= 19530 ~ str_glue("chr{Chromosome}:{Region}"),
as.numeric(bin_no) > 19530 & as.numeric(bin_no) <= 20547 ~ str_glue("chr{Chromosome}:{Region}"), as.numeric(bin_no) > 20547 & as.numeric(bin_no) <= 21512 ~ str_glue("chr{Chromosome}:{Region}"),
as.numeric(bin_no) > 21512 & as.numeric(bin_no) <= 22613 ~ str_glue("chr{Chromosome}:{Region}"), as.numeric(bin_no) > 22613 & as.numeric(bin_no) <= 23528 ~ str_glue("chr{Chromosome}:{Region}"),
as.numeric(bin_no) > 23528 & as.numeric(bin_no) <= 24307 ~ str_glue("chr{Chromosome}:{Region}"), as.numeric(bin_no) > 24307 & as.numeric(bin_no) <= 25228 ~ str_glue("chr{Chromosome}:{Region}"),
as.numeric(bin_no) > 25228 & as.numeric(bin_no) <= 26068 ~ str_glue("chr{Chromosome}:{Region}"), as.numeric(bin_no) > 26068 & as.numeric(bin_no) <= 26815 ~ str_glue("chr{Chromosome}:{Region}")))

#sample_wise_df_2$Chromosome <- NULL 
#sample_wise_df_2$Region <- NULL 
sample_wise_df_2$bin_no <- NULL

saveRDS(df_to_hmap_2, file = "results/k_50_table/k_50_filtered_region_mapped.rds")
readRDS("results/k_50_table/k_50_filtered_region_mapped.rds")

write.csv(df_to_hmap_2,"results/k_50_table/k_50_filtered_region_mapped.csv")
read.csv("results/k_50_table/k_50_filtered_region_mapped.csv")

write_feather(df_to_hmap_2,"results/k_50_table/k_50_filtered_region_mapped.feather")
write_feather(sample_wise_df_2,"results/k_50_table/k_50_sample_wise_df_region_mapped.feather")

write.csv(sample_wise_df_2,"results/k_50_table/k_50_sample_wise_df_region_mapped.csv")
read.csv("results/k_50_table/k_50_sample_wise_df_region_mapped.csv")