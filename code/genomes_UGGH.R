# Table 1 most common genera in human gut based on UGGH database

library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)

theme_set(theme_bw(base_size = 14))

genomes <- read.table('data/genomes-all_metadata.tsv', sep = '\t', header = T) %>% 
  dplyr::select(Genome, Genome_type, Completeness, Lineage, Sample_accession, Continent) %>% 
  dplyr::filter(Completeness > 90) %>%  
  tidyr::separate(Lineage, into=c('Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species'),
           sep=";") %>%  
  dplyr::mutate(Domain = str_remove_all(Domain, 'd__'), 
         Phylum = str_remove_all(Phylum, 'p__'), 
         Class = str_remove_all(Class, 'c__'), 
         Order = str_remove_all(Order, 'o__'), 
         Family = str_remove_all(Family, 'f__'), 
         Genus = str_remove_all(Genus, 'g__'), 
         Species = str_remove_all(Species, 's__'))

genomes

# By sample 
genus_n <- genomes %>% filter(Domain == 'Bacteria') %>% 
  group_by(Genus, Phylum, Sample_accession) %>% 
  reframe(ndist = n_distinct(Genome)) %>%  
  group_by(Phylum, Genus) %>% 
  reframe(n_genus_sample = mean(ndist))

genus_n %>% slice_max(n_genus_sample, n = 20) %>% 
  ggplot(aes(y = reorder(Genus, n_genus_sample), x = n_genus_sample, fill = Phylum )) +
  geom_col() +
  labs(x = '# genomes per sample', y = 'Genus')
ggsave('out/genomes_per_sample_UGGH.png')

# In the database 
genus_n2 <- genomes %>% filter(Domain == 'Bacteria') %>% 
  group_by(Genus, Phylum) %>% 
  reframe(ndist = n_distinct(Genome)) 

genus_n2 %>% slice_max(ndist, n = 20) %>% 
  ggplot(aes(y = reorder(Genus, ndist), x = ndist, fill = Phylum)) +
  geom_col() +
  labs(x = '# genomes in database', y = 'Genus')
ggsave('out/genomes_per_database_UGGH.png')

# Species 
species <- genomes %>% filter(Domain == 'Bacteria', Continent == 'Europe') %>% 
  group_by(Species, Phylum) %>% 
  reframe(ndist = n_distinct(Genome)) 
  # group_by(Phylum, Species) %>% 
  # reframe(n_genus_sample = mean(ndist))

species %>% slice_max(ndist, n = 20) %>% 
  ggplot(aes(y = reorder(Species, ndist), x = ndist, fill = Phylum )) +
  geom_col() +
  labs(x = '# genomes per database', y = 'Species')

species %>% slice_max(n_genus_sample, n = 20) %>% 
  ggplot(aes(y = reorder(Species, n_genus_sample), x = n_genus_sample, fill = Phylum )) +
  geom_col() +
  labs(x = '# genomes per sample', y = 'Species')


