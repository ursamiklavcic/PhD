# Sporulation ability determination 

# Library
library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)
library(stringr)
library(tibble)
library(purrr)

set.seed(96)
theme_set(theme_bw())

# Import results from code/sporulation_ability/03_sporeability.sh

blast_results <- read.table('data/sporulation_ability/mpa_blast_results.tsv', sep = '\t', header = TRUE) %>% 
  mutate(genome_id = substr(genome_id, 1, 15)) %>% 
  rename('locus_tag' = 'gene_name')

# Sporulation genes info (weight etc)
gene_info <- read.table('data/sporulation_ability/gene_info.csv', sep = ';', header = TRUE) %>% 
  select(locus_tag, gene_name, weight)

# Import results from 02_find_genomes_from_species.sh
genome_info <- read.table('data/sporulation_ability/mpa_accession_manual.tsv', sep = '\t', header = TRUE)


abund <- read_tsv('~/projects/longitudinal_shotgun/data/metaphlan_abundance_table.txt', comment = '#') %>%
  rename_with(~ str_remove(., '^profiled_'), starts_with('profiled_')) %>% 
  pivot_longer(-clade_name) %>% 
  filter(grepl('s__', clade_name), !grepl('t__', clade_name)) 
  # mutate(clade_name = str_remove_all(clade_name, '[a-zA-Z]__')) %>%
  # separate(clade_name, into=c('Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species', 'SGB'),
  #          sep="\\|") %>% 
  # mutate(Phylum = ifelse(Phylum == 'p__Firmicutes', 'p__Bacillota', Phylum), 
  #        Domain = str_remove_all(Domain, 'k__'), 
  #        Phylum = str_remove_all(Phylum, 'p__'), 
  #        Class = str_remove_all(Class, 'c__'), 
  #        Order = str_remove_all(Order, 'o__'), 
  #        Family = str_remove_all(Family, 'f__'), 
  #        Genus = str_remove_all(Genus, 'g__'), 
  #        Species = str_remove_all(Species, 's__'), 
  #        SGB = str_remove_all(SGB, 't__')) 
  

# Sporulation ability based on Browne et al. 2017
blastpre <- blast_results %>% 
  left_join(genome_info, by = 'genome_id', relationship = 'many-to-many') %>% 
  left_join(gene_info, by = 'locus_tag', relationship = 'many-to-many') %>% 
  filter(evalue < 10e-5 & identity > 30) 

blast <- blastpre %>% 
  group_by(genome_id, clade_name) %>% 
  reframe(spore_genes = n_distinct(locus_tag), 
          raw_score = sum(weight)) %>%  
  mutate(sporulation_score = raw_score/(max(raw_score)), 
         spore_former = ifelse(sporulation_score >= 0.5, TRUE, FALSE))

# Plot sporulation signature 
ggplot(blast, aes(x = spore_genes, y = sporulation_score)) +
  geom_point(size = 3) +
  geom_vline(xintercept = 33) +
  geom_hline(yintercept = 0.5)
ggsave('out/sporulation/sporulation_score_n_genes_.png')


# 
# Sporulation ability based on Browne et al. 2021
gene_counts <- blastpre %>%  
  group_by(genome_id, clade_name) %>% 
  reframe(n_gene = n_distinct(locus_tag)) 

# Family distributions 
gene_counts %>% 
  separate(clade_name, into=c('Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species'), sep="\\|") %>% 
  ggplot(aes(x = n_gene)) +
  geom_histogram() +
  geom_vline(xintercept = 33) +
  facet_wrap(~Family, scales = 'free_y')
ggsave('out/sporulation/spore_genes_family.png', width = 29, height = 15, units = 'cm')

spore_ability <- blastpre %>% 
  separate(clade_name, into=c('Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species'), sep="\\|") %>% 
  group_by(genome_id) %>% 
  reframe(PA = any(gene_name == 'spo0A'), 
          n_genes = n_distinct(locus_tag)) %>% 
  mutate(sporulation_ability = ifelse(PA == TRUE & n_genes  >= 33,  "Spore-former",  "Non-spore-former"))

spore_ability %>%  group_by(PA) %>%  
  reframe(n = n_distinct(genome_id))

spore_ability %>% 
  group_by(sporulation_ability) %>% 
  reframe(n = n_distinct(genome_id))

spore_ability
# 
# Save
spores <- spore_ability %>% 
  left_join(genome_info, by = c('genome_id')) 

write.table(spores, file='data/sporulation_ability/sporulation_ability2021.tsv',
            quote=FALSE, row.names = FALSE, sep='\t')


