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

# Sporulation ability 
# Which species have I not yet defined as sporeforming?? 
sporeforming <- read.table('../spore_evo_dynamics/data/sporulation_signature/sporulation_ability_all_withTAX.tsv', sep = '\t', header = TRUE)
colnames(sporeforming) <- c('genome_id', 'species_id', 'Phylum', 'Order', 'Class', 'Family', 'Genus', 'Species', 'n_spore_genes', 'sporeforming')

sporeforming <- mutate(sporeforming, Phylum = str_replace(Phylum, "^Firmicutes", "Bacillota"))

# GTDB. Get the species
gtdb <- read.table('~/projects/longitudinal_shotgun/data/gtdb_merged.txt', sep = '\t', header = TRUE) %>% 
  separate(clade_name, into=c('Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species'), sep=";") %>% 
  #filter(!is.na(Species)) %>% 
  pivot_longer(-c('Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species')) %>% 
  mutate(Domain = str_remove_all(Domain, 'd__'), 
         Phylum = str_remove_all(Phylum, 'p__'), 
         Class = str_remove_all(Class, 'c__'), 
         Order = str_remove_all(Order, 'o__'), 
         Family = str_remove_all(Family, 'f__'), 
         Genus = str_remove_all(Genus, 'g__'), 
         Species = str_remove_all(Species, 's__'), 
         name = str_remove_all(name, 'gtdbdd_'), 
         name = str_remove_all(name, '.txt')) %>% 
  filter(Domain == 'Bacteria', !is.na(Phylum), !is.na(Class), 
         !is.na(Order), !is.na(Family), !is.na(Genus), !is.na(Species)) %>% 
  mutate(PA = ifelse(value > 0, 1, 0)) %>% 
  left_join(sporeforming, by = c('Phylum', 'Order', 'Class', 'Family', 'Genus', 'Species'))

# We have defined sporeforming ability for the ones found with MIDAS not metaphlan. Now I want to add those that are missing! 
species_to_define <- filter(gtdb, is.na(sporeforming)) %>% 
  select(Species) %>% 
  unique()

write_csv(species_to_define, 'species_to_define.csv')

# Downloaded GTDB metadata 
# go to advanced search and select all GTDB representative genomes (download and transfer to HPC) 
gtdb_meta <- read.table('gtdb-adv-search.tsv', sep = '\t', header = TRUE, fill = TRUE, quote="") %>% 
  separate(GTDB.Taxonomy, into=c('Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species'), sep=";") %>%
  mutate(Species = str_remove_all(Species, 's__')) %>% 
  select(Accession, Species)

species_to_define_GCF <- left_join(species_to_define, gtdb_meta, by = 'Species')
write_tsv(species_to_define_GCF, 'species_to_define_GCF.tsv')

# sporeability.sh :  
# - Download genomes 
# - Find sporulation genes tBLASTn
# - Import and count the number of sporulation genes! 

blast_results <- 
  # First from MIDAS work
  read_tsv('~/projects/spore_evo_dynamics/data/sporulation_signature/all_blast_results.tsv') %>% 
  rename('locus_tag' = 'gene_name') %>% 
  mutate(genome_id = sub("^(([^_]+)_([^_]+))_.*", "\\1", genome_id)) %>% 
  # missing ones from now
  rbind(read_tsv('~/projects/spore_evo_dynamics/data/sporulation_signature/additional_mpa_blast_results.tsv') %>% 
          rename('locus_tag' = 'gene_name') %>% 
          mutate(genome_id = sub("^(([^_]+)_([^_]+))_.*", "\\1", genome_id)))

# Sporulation genes info (weight etc)
gene_info <- read.table('longitudinal_shotgun/data/gene_info.csv', sep = ';', header = TRUE) 

# Genomes from MIDAS 
g1  <- read_tsv('~/projects/spore_evo_dynamics/data/sporulation_signature/species_for_analysis.tsv', col_names = FALSE) %>%
  select(X1, X6) %>% 
  rename('genome_id' = 'X1', 'tax' = 'X6') 

# GTDB representative genomes
g2 <- read.table('longitudinal_shotgun/data/gtdb-adv-search.tsv', sep = '\t', header= TRUE, fill = TRUE, quote = '') %>%  
  select(Accession, GTDB.Taxonomy) %>% 
  rename('genome_id' = 'Accession', 'tax' = 'GTDB.Taxonomy')

# All genomes from GTDB
g3 <- read.table('longitudinal_shotgun/data/bac120_taxonomy_r220.tsv.gz', sep = '\t', header = FALSE) %>%  
  rename('genome_id' = 'V1', 'tax' = 'V2') %>% 
  mutate(genome_id = str_remove_all(genome_id, 'RS_|GB_')) 
  # separate(V2, into=c('Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species'), sep=";") 

# Manually found genomes! 
# genome_info <- rbind(g1, g2, g3) %>% unique()
# blast <- blast_results %>% 
#   left_join(genome_info, by = 'genome_id', relationship = 'many-to-many')
# 
# # The ones that are still not defined. Search manually in GTDB !
# g <- filter(blast, is.na(Domain)) %>% 
#   select(genome_id) %>% 
#   unique() %>%  
#   mutate(n = n_distinct(genome_id))
# write_csv(g, 'g.csv')  
  
g4 <- read_csv2('longitudinal_shotgun/data/genomes.csv') 

# All genomes taxonomy 
genome_info <- rbind(g1, g2, g3, g4) %>% distinct(genome_id, .keep_all = TRUE) 
  

# Sporulation ability based on Browne et al. 2017
blastpre <- blast_results %>% 
  left_join(genome_info, by = 'genome_id', relationship = 'many-to-many') %>% 
  left_join(select(gene_info, locus_tag, gene_name, weight), by = 'locus_tag', relationship = 'many-to-many') %>% 
  filter(evalue < 10e-5 & identity > 30) 

# Check if all genomes have taxonomy classification! 
blastpre %>%  filter(is.na(tax))

blast <- blastpre %>% 
  group_by(genome_id, tax) %>% 
  reframe(spore_genes = n_distinct(locus_tag), 
          raw_score = sum(weight)) %>%  
  mutate(sporulation_score = raw_score/(max(raw_score)), 
         spore_former = ifelse(sporulation_score >= 0.5, TRUE, FALSE))

# Save
# tsv <- blast %>%
#   mutate(sporeforming = ifelse(scaled_score > 0.499, TRUE, FALSE)) %>%
#   select(genome_id, sporeforming)
# 
# write.table(tsv, file='longitudinal_shotgun/data/sporulation_ability2017.tsv',
#             quote=FALSE, row.names = FALSE, sep='\t')

# Plot sporulation signature 
ggplot(blast, aes(x = spore_genes, y = sporulation_score)) +
  geom_point(size = 3) +
  geom_vline(xintercept = 33) +
  geom_hline(yintercept = 0.5)
ggsave('out/sporulation/sporulation_score_n_genes_.png')


# 
# Sporulation ability based on Browne et al. 2021
gene_counts <- blastpre %>%  
  group_by(genome_id, tax) %>% 
  reframe(n_gene = n_distinct(gene_name)) 

# Family distributions 
gene_counts %>% 
  separate(tax, into=c('Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species'), sep=";") %>% 
  ggplot(aes(x = n_gene)) +
  geom_histogram() +
  geom_vline(xintercept = 33) +
  facet_wrap(~Family, scales = 'free_y')
ggsave('out/sporulation/spore_genes_family.png')

spo0A <- blastpre %>% 
  separate(tax, into=c('Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species'), sep=";") %>% 
  group_by(genome_id, Family) %>% 
  reframe(PA = any(gene_name == 'spo0A'))

spo0A %>%  group_by(PA) %>%  
  reframe(n = n_distinct(genome_id))

family_summary <- gene_counts %>%
  separate(tax, into=c('Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species'), sep=";") %>% 
  group_by(Family) %>%
  reframe(n_genomes = n_distinct(genome_id), 
          n_genes = sum(n_gene)/n_genomes) 

spore_ability <- spo0A %>% 
  left_join(family_summary, by = 'Family', relationship = 'many-to-many') %>% 
  mutate(sporulation_ability = ifelse(PA == TRUE & n_genes  >= 33,  "Spore-former",  "Non-spore-former"))

spore_ability %>% 
  group_by(sporulation_ability) %>% 
  reframe(n = n_distinct(genome_id))

spore_ability
# 
# Save
spores <- select(spore_ability, genome_id, n_genes, sporulation_ability) %>% 
  left_join(genome_info, by = c('genome_id')) 

write.table(spores, file='longitudinal_shotgun/data/sporulation_ability2021.tsv',
            quote=FALSE, row.names = FALSE, sep='\t')


