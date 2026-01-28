# Quntify genes in metagenomic samples 

library(tidyr)
library(dplyr)
library(ggplot2)
library(stringr)

theme_set(theme_bw(base_size = 14))

# Find genes specific for each strain
# We will use this genes to map metagenomic reads back to them and get 
# a better representation of which strains were present when within individual H at each time point

genes <- read.csv('~/projects/C_saudiense/004_prokka/roary/gene_presence_absence.csv') %>% 
  pivot_longer(names_to = 'samples', values_to = 'PA', cols = starts_with('H00')) %>%  
  mutate(PA_01 = ifelse(PA != '', 1, 0), 
         samples = str_remove_all(samples, ".fasta"), 
         samples = str_remove_all(samples, "_hybrid"))

metadata_genomes <- read.table('~/projects/C_saudiense/metadata_genomes.tsv', sep = '\t', header = T)
genes <- left_join(genes, metadata_genomes, by = 'samples')

no_strains <- group_by(genes, STRAINS) %>%  
  reframe(n = n_distinct(samples))
no_strains


strain_genes <- genes %>%  
  filter(PA_01 == 1) %>% 
  group_by(Gene, Annotation) %>%  
  reframe(n = n_distinct(STRAINS), 
          who = list(unique(STRAINS)), 
          ID = PA) %>% 
  filter(n== 1) %>%  
  mutate(STRAIN = unlist(who)) %>%  
  select(-c(who, n))

# write.csv2(strain_genes, 'out/specific_genes_per_strain.csv')

# For thesis
thesis_strain <- strain_genes %>%  group_by(STRAIN) %>% 
  reframe(Gene = list(unique(Gene))) %>% 
  unnest(cols = c(Gene)) %>% 
  left_join(select(strain_genes, Gene, Annotation) %>% distinct(), by = 'Gene')

write.csv2(thesis_strain, 'out/thesis_strain_genes.csv')

# Number unique genes per strain 
strain_genes %>% #filter(Annotation == 'hypothetical protein') %>% 
  group_by(STRAIN) %>% 
  reframe(n = n_distinct(Gene))


# RUN LARA's SCRIPTS TO OBTAIN GENE SEQUENCES, later this code 

# value is already normalized by the gene length! 
genes_tpm <- read.table('~/projects/C_saudiense/006_quantify_genes_Lara/combine_tpm.txt', header = T) %>% 
  pivot_longer(-c(Gene, Gene_ID))

metadata <- read.table('~/projects/longitudinal_shotgun/data/metadata.csv', sep = ';', header = TRUE)

strain_genes <- read.table('~/projects/C_saudiense/out/specific_genes_per_strain.csv', sep = ';', header = T) 

genes_meta <- left_join(genes_tpm, select(strain_genes, ID, STRAIN), by = join_by('Gene_ID' == 'ID')) %>% 
  left_join(metadata, by = join_by('name' == 'Group')) %>% 
  # C. saudiense was detected in sample MA003 ONLY! 
  filter(name != 'MA003')

# Strains in person H or A ? 
# 
genes_meta %>%  
  filter(value > 0) %>% 
  group_by(STRAIN, person) %>%  
  reframe(n = n_distinct(Gene))

# Genes detected in both participants? 
both <- genes_meta %>%  
  filter(value > 0) %>%  
  group_by(Gene, person, STRAIN) %>%  
  reframe(mean_value = mean(value, na.rm = T))
  
both2 <- filter(both, person == 'H') %>% 
  full_join(filter(both, person == 'A'), by = c('Gene', 'STRAIN')) %>% 
  mutate(true = ifelse(mean_value.x < mean_value.y, TRUE, FALSE))

both2 %>% filter(!is.na(mean_value.x), !is.na(mean_value.y)) %>% 
  group_by(STRAIN) %>% 
  reframe(n = n_distinct(Gene))


# Are genes truly not detected in samples of participant A and only in samples of participant H? 
genes_meta %>%  
  filter(value > 0) %>% 
  group_by(Gene, person) %>% 
  reframe(n = n_distinct(name)) %>%  
  ggplot(aes(x = Gene, y = n, color = person)) +
  geom_point(size = 3) +
  theme(axis.text.x = element_text(angle = 90))

genes_specific <- genes_meta %>%  
  filter(value > 0) %>% 
  group_by(Gene) %>% 
  reframe(n = n_distinct(person)) %>%  
  group_by(n) %>%  
  reframe(n_genes = n_distinct(Gene))
genes_specific
# 1 = 1113
# 2 = 85

# What are these genes that can be found in participant A, even though C. saudiese is not? 
genes_shared <- genes_meta %>%  
  filter(value > 0) %>% 
  group_by(Gene, Gene_ID) %>% 
  reframe(n = n_distinct(person), 
          who_person = list(unique(person))) %>%  
  group_by(n) %>%  
  reframe(n_genes = n_distinct(Gene), 
          who = list(unique(Gene_ID))) %>% 
  filter(n == 2) %>% 
  pull(unlist(who))

genes_shared_list <- filter(genes, PA %in% genes_shared[[1]])

# There is one gene only detected within person A which is it?
genes_meta %>%  
  filter(value > 0) %>% 
  group_by(Gene, Gene_ID) %>% 
  reframe(n = n_distinct(person), 
          who_person = list(unique(person))) %>% 
  filter(n == 1) %>% 
  mutate(person = unlist(who_person)) %>%  
  filter(person == 'A')
# BOKAGMOO_01020     

# Do I filter our the genes that are not specific? YES
genes_filt <- genes_meta %>%  
  filter(value > 0, !Gene_ID %in% genes_shared[[1]], Gene_ID != 'BOKAGMOO_01020') 

genes_filt %>% 
  ggplot(aes(x = day, y = value, color = person)) +
  geom_point(size = 4) +
  facet_wrap(person~biota, scales = 'free') +
  labs(x = 'Day', y = 'TPM')

# Number of genes we started with
no_genes <- genes_meta %>% group_by(STRAIN) %>% 
  reframe(n_genes = n_distinct(Gene))

# How many genes were found per strain and sample?
genes_filt %>%  
  group_by(name, STRAIN) %>% 
  reframe(n_genes_found = n_distinct(Gene)) %>% 
  left_join(metadata, by = join_by('name' == 'Group')) %>% 
  ggplot(aes(x = as.factor(day), y = n_genes_found, fill = STRAIN)) +
  geom_col() +
  facet_wrap(person ~ biota, scales = 'free')

# As different STRAINS have different number of genes, I have to normalize TPM (value) 
# to take into account the number of genes per strain 
no_genes2 <- genes_filt %>% 
  group_by(STRAIN) %>%  
  reframe(n_genes = n_distinct(Gene)) 

events <- read.table('~/projects/thesis/data/extreme_event_data.csv', sep = ',', header =T) %>% 
  filter(person == 'H')

genes_filt %>% 
  left_join(no_genes2, by = 'STRAIN') %>% 
  group_by(name, STRAIN) %>% 
  reframe(norm_TPM = sum(value)/n_genes) %>% 
  filter(norm_TPM > 0) %>% 
  left_join(metadata, by = join_by('name' == 'Group')) %>% 
  mutate(biota = ifelse(biota == 'bulk microbiota', 'untreated sample', biota)) %>% 
  ggplot(aes(x = day, y = norm_TPM, color = STRAIN)) +
  geom_rect(data = events, aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = Inf, fill = extremevent_type), inherit.aes = FALSE,
            alpha = 0.6) +
  scale_fill_manual(values = c('white','#d94343', '#d98e43', '#f1f011', '#0c9910', '#3472b7', '#7934b7', '#b73485', '#0f5618')) +
  geom_line(linewidth = 2) +
  geom_point(size = 3) +
  facet_wrap(~ biota, scales = 'free') +
  labs(x = 'Day', y = '# reads per million reads normalized by\nlength of gene and number of genes per strain', 
       fill = 'Event', color = 'Strain')
ggsave('out/quantify_genes_H.svg')

# log scale 
genes_filt %>% 
  left_join(no_genes2, by = 'STRAIN') %>% 
  group_by(name, STRAIN) %>% 
  reframe(norm_TPM = sum(value)/n_genes) %>% 
  filter(norm_TPM > 0) %>% 
  left_join(metadata, by = join_by('name' == 'Group')) %>% 
  mutate(biota = ifelse(biota == 'bulk microbiota', 'untreated sample', biota)) %>% 
  ggplot(aes(x = day, y = norm_TPM, color = STRAIN)) +
  geom_rect(data = events, aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = Inf, fill = extremevent_type), inherit.aes = FALSE,
            alpha = 0.6) +
  scale_fill_manual(values = c('white','#d94343', '#d98e43', '#f1f011', '#0c9910', '#3472b7', '#7934b7', '#b73485', '#0f5618')) +
  geom_line(linewidth = 2) +
  geom_point(size = 3) +
  scale_y_log10() +
  facet_wrap(~ biota, scales = 'free') +
  labs(x = 'Day', y = '# reads per million reads normalized by\nlength of gene and number of genes per strain [log10]', 
       fill = 'Event', color = 'Strain')
ggsave('out/quantify_genes_H_log10.svg')


# This is used in thesis!!!
# What if instead of working with TPM, I work with count and normalize it by length of gene and number of genes per strain? 
genes_count_pre <- read.table('~/projects/C_saudiense/006_quantify_genes_Lara/combined_counts.txt', sep = '', header = T) %>%  
  pivot_longer(-Gene_ID)

genes_norm_count <- left_join(genes_count_pre, select(strain_genes, ID, STRAIN), by = join_by('Gene_ID' == 'ID')) %>%
  # Filter genes found in more than just individual H
  filter(!Gene_ID %in% genes_shared[[1]], Gene_ID != 'BOKAGMOO_01020') %>% 
  # join info on genes length
  left_join(select(genes, PA, Avg.group.size.nuc), by = join_by('Gene_ID' == 'PA')) %>%
  # normalize by gene length
  mutate(norm_count = value*150/Avg.group.size.nuc) %>% 
  left_join(no_strains, by = join_by('STRAIN' == 'STRAINS')) %>% 
  group_by(name, STRAIN) %>% 
  # normalize by number of genes in a strain
  reframe(sum_norm_count = sum(norm_count)/n) %>% 
  left_join(metadata, by = join_by('name' == 'Group')) %>% 
  filter(name != 'MA003', !is.na(sum_norm_count))


genes_norm_count %>% 
  filter(person == 'H') %>% 
  mutate(biota = ifelse(biota == 'bulk microbiota', 'untreated sample', biota)) %>% 
  ggplot(aes(x = day, y = sum_norm_count, color = STRAIN)) +
  geom_rect(data = events, aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = Inf, fill = extremevent_type), inherit.aes = FALSE,
            alpha = 0.6) +
  scale_fill_manual(values = c('white','#d94343', '#d98e43', '#f1f011', '#0c9910', '#3472b7', '#7934b7', '#b73485', '#0f5618')) +
  geom_line(linewidth = 2) +
  geom_point(size = 3) +
  facet_wrap( ~ biota, scales = 'free') +
  #scale_y_log10() +
  labs(x = 'Day', y = '# reads normalized by the gene length\nand number of genes per strain', fill = 'Event')
ggsave('out/counts_normalized_by_gene_length.svg', dpi = 600)

# log 10
genes_norm_count %>% 
  filter(person == 'H') %>% 
  mutate(biota = ifelse(biota == 'bulk microbiota', 'untreated sample', biota)) %>% 
  ggplot(aes(x = day, y = sum_norm_count, color = STRAIN)) +
  geom_rect(data = events, aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = Inf, fill = extremevent_type), inherit.aes = FALSE,
            alpha = 0.6) +
  scale_fill_manual(values = c('white','#d94343', '#d98e43', '#f1f011', '#0c9910', '#3472b7', '#7934b7', '#b73485', '#0f5618')) +
  geom_line(linewidth = 2) +
  geom_point(size = 3) +
  facet_wrap( ~ biota, scales = 'free') +
  scale_y_log10() +
  labs(x = 'Day', y = '# reads normalized by the gene length\nand number of genes per strain [log10]', fill = 'Event')
ggsave('out/counts_normalized_by_gene_length_log10.png', dpi = 600)
ggsave('out/counts_normalized_by_gene_length_log10.svg', dpi = 600)

