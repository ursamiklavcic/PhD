# Sporulation genes in metagenomes 

library(tidyr)
library(stringr)
library(dplyr)
library(ggplot2)

metadata <- read.table('data/metadata.csv', sep = ';', header = T)
metadata_genes <- read.table('data/sporulation_ability/gene_info.csv', sep = ';', header = T)
gene_length <- read.table('data/sporulation_ability/metagenomes/index/gene_lengths.txt', header = F)

gene_metagenomes <- read.table('data/sporulation_ability/metagenomes/data/combined_counts.txt', header = T) %>% 
  left_join(metadata_genes, by = 'Gene_ID') %>% 
  left_join(gene_length, by = join_by('Gene_ID' == 'V1')) %>% 
  pivot_longer(names_to = 'Group', values_to = 'count', cols = starts_with(c('M', 'S'))) %>% 
  filter(count > 0) %>% 
  # Normalize the number of counts by gene length 
  mutate(norm_count = (count*150)/V2) %>% 
  select(Group, locus_tag, norm_count) %>% 
  distinct() %>% 
  left_join(metadata, by = 'Group')

# Number of sporulation genes in samples of bulk microbiota and sporobiota? 
gene_metagenomes %>% 
  ggplot(aes(x = Group, y = n_distinct(locus_tag), fill = biota)) +
  geom_col()
