# PStrain with sporulation ability 






# some old stuff - look thhrough 



species_rel <- read_tsv('data/species/species_relative_abundance.tsv.gz') %>% 
  pivot_longer(names_to = 'sample_name', values_to = 'rel_abund', cols = -species_id)

blast_rel <- blast %>%  
  full_join(species_rel, by = 'species_id')
metadata <- read_csv2('data/filtered_metadata.csv')

blast_rel %>% 
  filter(spore_genes > 0 & rel_abund > 0) %>% 
  ggplot(aes(x = spore_genes, y = rel_abund, color = phylum)) +
  geom_jitter(size = 2, alpha = .8) +
  geom_vline(xintercept = 33) +
  scale_y_log10() +
  scale_color_manual(values = c('#d94343', '#d98e43', '#f1f011', '#0f5618','#0c9910',
                                '#23bf6c', '#3472b7', '#7934b7', '#b73485')) +
  labs(x = '# spore genes', y = 'log10(relative abundance', color = '') 
ggsave('plots/no_spore_genes_relabund.tiff', dpi = 600)

blast_prev <- blast_rel %>%  
  mutate(PA = ifelse(rel_abund > 0, 1, 0)) %>% 
  left_join(metadata, by = 'sample_name') %>% 
  group_by(person, biota, genome_id, spore_genes) %>% 
  reframe(prevalence = (sum(PA)/n())*100) %>% 
  filter(!is.na(genome_id))

max(blast_prev$prevalence)

blast_prev %>%
  filter(biota == 'bulk microbiota') %>% 
  ggplot(aes(x = spore_genes, y = prevalence, color = person)) +
  geom_jitter() +
  geom_vline(xintercept = 33) +
  labs(x = '# spore genes', y = 'prevalence (%)', 
       caption = 'Prevalence of each species within each individual thriugh time')
ggsave('plots/no_spore_genes_prevalence.tiff', dpi=600)

blast_rel %>%  
  left_join(metadata, by = 'sample_name') %>% 
  filter(rel_abund > 0) %>% 
  ggplot(aes(x = spore_genes, y = rel_abund, color = biota)) +
  geom_point()

# Count the number of species in the ethanol treated samples have more than 33 sporulation genes
# vs how many in the bulk microbiota 

blast_rel %>%  
  left_join(metadata, by = 'sample_name') %>%
  mutate(sporeforming = ifelse(spore_genes > 32, TRUE, FALSE)) %>%
  filter(rel_abund > 0) %>% 
  group_by(person, biota, sporeforming) %>%  
  reframe(n = n_distinct(species)) %>%  
  filter(!is.na(sporeforming) & sporeforming == FALSE) %>%  
  pivot_wider(names_from = biota, values_from = n)

### 
# For all species present in the dataset 
# Sporulation ability of all species detected in our dataset 
blast_results_all <- read_tsv('data/sporulation_signature/all_blast_results.tsv') %>% 
  rename('locus_tag' = 'gene_name') %>% 
  mutate(genome_id = sub("^(([^_]+)_([^_]+))_.*", "\\1", genome_id))

genome_info <- read_tsv('data/genomes_metadata.tsv', col_names = FALSE) %>% 
  rename('genome' = 'X1', 'species_id' = 'X2', 
         'genome_id' = 'X3', 'genome_is_representative' = 'X4', 
         'species' = 'X5', 'tax' = 'X6', 
         'ANI_circumscription_radius' = 'X7', 
         'Mean_intraspecies_ANI' = 'X8', 
         'Min_intraspecies_ANI' = 'X9', 
         'Mean_intraspecies_AF' = 'X10', 
         'Min_intraspecies_AF' = 'X11', 
         'No_clustered_genomes' = 'X12') %>% 
  mutate(tax = str_remove_all(tax, '[a-zA-Z]__'), 
         species = str_remove_all(species, '[s]__'))

blast_all <- blast_results_all %>% 
  left_join(select(genome_info, genome_id, species, tax, species_id), by = 'genome_id') %>% 
  separate(tax, into=c('domain', 'phylum', 'class', 'order', 'family', 'genus'), sep=";") %>% 
  left_join(select(gene_info, locus_tag, gene_name, description), by = 'locus_tag', relationship = 'many-to-many') %>% 
  filter(evalue < 10e-5 & identity > 30) %>% 
  group_by(genome_id, species_id, phylum, order, class, family, genus, species) %>% 
  reframe(spore_genes = n_distinct(locus_tag)) 

blast_all %>% 
  ggplot(aes(x = reorder(species, spore_genes), y = spore_genes, )) +
  geom_point(size=1) +
  scale_y_continuous(breaks = c(0,5,10, 15,20,25,30,33,35,40,50,60,66)) +
  theme(axis.text.x = element_blank(),
        panel.grid.major.x = element_blank() ,
        # explicitly set the horizontal lines (or they will disappear too)
        panel.grid.major.y = element_line( size=.1, color="black" )) +
  labs(x = 'Species reorderd by # spore genes present', y = '# spore genes')
ggsave('plots/no_sporulation_genes_per_species_all.tiff', dpi=600)

# Histogram 
blast_all %>% 
  ggplot(aes(x = spore_genes)) +
  geom_bar() +
  scale_x_continuous(breaks = c(0,5,10, 15,20,25,30,33,35,40,50,60,66)) +
  geom_vline(xintercept = 33) +
  labs(x = '# spore genes', y = '# species')
ggsave('plots/no_sporulation_genes_per_species_all_histo.tiff', dpi=600)


tsv <- blast_all %>% 
  mutate(sporeforming = ifelse(spore_genes > 33, TRUE, FALSE)) %>% 
  select(genome_id, species_id, sporeforming)

write.table(tsv, file='data/sporulation_signature/sporulation_ability_all.tsv', 
            quote=FALSE, row.names = FALSE, sep='\t')



# Are sporeforming bacteria found in ethanol treated samples? 
ethanol_resistancy <- read_tsv('data/ethanol_resist.tsv') 

blast_rel %>%  
  left_join(ethanol_resistancy, by = 'species_id') %>% 
  ggplot(aes(x = spore_genes, y = reorder(species, spore_genes), color = ethanol_resistant)) +
  geom_point(size = 3) +
  facet_wrap(~ethanol_resistant)

blast_rel %>%  
  left_join(ethanol_resistancy, by = 'species_id') %>% 
  ggplot(aes(x = rel_abund, y = reorder(species, rel_abund), color = ethanol_resistant)) +
  geom_point(size = 3) +
  scale_x_log10()

