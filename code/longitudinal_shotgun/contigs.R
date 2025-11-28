library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)
library(stringr)
library(lubridate)
library(ggpubr)

set.seed(96)
theme_set(theme_bw(base_size = 14))

# Metadata 
metadata <- read.table('data/metadata.csv', header= TRUE, sep = ';') %>%
  mutate(date = dmy(date)) 

# Analysis can be run on contigs, ORF or bins 
# Contigs 
contigs_pre <- read.table('~/projects/longitudinal_shotgun/data/sqm_tables/19.sqmMicrobiota.contigtable', sep = '\t', header = TRUE, comment.char = '#') %>%
  # remove contigs less than 1000 bp long and whole disparity is more than 10%
  filter(Length > 1000 & Disparity < 0.1)

# # ORFs
# orfs <- read_csv('data/tables/13.sqmMicrobiota.orftable')
# 
# # Bins
# bins <- read_csv('data/tables/18.sqmMicrobiota.bintable')

# Taxonomic composition based on contigs, taxonomy of which has been derived by the LCA algorithm from ORFs on this contig! 
contigs <- contigs_pre %>%
  # rename some columns that have spaces in names
  rename(contigID = `Contig ID`, NoGenes =  `Num genes`) %>%
  # # select onyl columns that I will need
  # select(contigID, Tax, Disparity, Length, NoGenes, starts_with('TPM')) %>%
  # # rename so that i have sample names
  # rename_with(~ str_remove(., 'TPM '), starts_with('TPM')) %>%
  # remove prefix and separate Tax column
  mutate(Tax = str_replace_all(Tax, '[a-zA-Z]_', '')) %>%
  separate(Tax, into=c('Kingdom', 'Clade', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species'),
           sep=";") %>%
  mutate(Clade = ifelse(is.na(Clade), Kingdom, Clade), 
         Phylum = ifelse(is.na(Phylum), Clade, Phylum), 
         Class = ifelse(is.na(Class), Phylum, Class), 
         Order = ifelse(is.na(Order), Class, Order), 
         Family = ifelse(is.na(Family), Order, Family), 
         Genus = ifelse(is.na(Genus), Family, Genus), 
         Species = ifelse(is.na(Species), Genus, Species))

saveRDS(contigs, 'data/intermediate/contigs.RDS')

contigs <- readRDS('~/projects/longitudinal_shotgun/data/intermediate/contigs.RDS')

kingdom_contigs <- contigs %>% 
  mutate(TPM = select(., `TPM MA001`:`TPM MI012`) %>% rowSums(na.rm = TRUE), x = 'X') %>%  
  group_by(Kingdom, x) %>% 
  reframe(TPM = sum(TPM)) %>%  
  mutate(rel_abund = TPM/sum(TPM)*100, 
         Kingdom = ifelse(is.na(Kingdom), 'UNCLASSIFIED', Kingdom))

conitgs_domain <- kingdom_contigs %>% 
  mutate(Kingdom = factor(Kingdom, levels = c('Archaea', 'Bacteria', 'Eukaryota', 'Viruses', 'UNCLASSIFIED')), 
         x = 'DIAMOND + LCA') %>% 
  ggplot(aes(x = rel_abund, y = x, fill = Kingdom)) +
  geom_col() +
  scale_fill_manual(values = c('#4E9E23', '#29BCE3', '#D13D21', '#FAE528', 'grey')) +
  labs(y = '', x = 'Relative abundance [%]', fill = 'Domain')+
  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank())



abund <- read_tsv('~/projects/longitudinal_shotgun/data/metaphlan_abundance_table.txt', comment = '#') %>%
  rename_with(~ str_remove(., '^profiled_'), starts_with('profiled_')) %>%
  #mutate(clade_name = str_remove_all(clade_name, '[a-zA-Z]__')) %>%
  separate(clade_name, into=c('Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species', 'SGB'),
           sep="\\|") %>% 
  mutate(Phylum = ifelse(Phylum == 'p__Firmicutes', 'p__Bacillota', Phylum), 
         Domain = str_remove_all(Domain, 'k__'), 
         Phylum = str_remove_all(Phylum, 'p__'), 
         Class = str_remove_all(Class, 'c__'), 
         Order = str_remove_all(Order, 'o__'), 
         Family = str_remove_all(Family, 'f__'), 
         Genus = str_remove_all(Genus, 'g__'), 
         Species = str_remove_all(Species, 's__'), 
         SGB = str_remove_all(SGB, 't__')) %>% 
  select(-MC013)

domain <- filter(abund, is.na(Phylum)) %>% 
  select(-c(Phylum, Class, Order, Family, Genus, Species, SGB)) %>% 
  pivot_longer(-Domain) %>%  
  left_join(metadata, by = join_by('name' == 'Group')) %>% 
  filter(biota == 'bulk microbiota') %>% 
  mutate(rel_abund = value /sum(value) * 100) %>% 
  group_by(Domain) %>% 
  reframe(rel_abund = sum(rel_abund)) %>%  
  mutate( x = 'MetaPhlAn') %>%  
  ggplot(aes(x = rel_abund, y = x, fill = Domain)) +
  geom_col() +
  scale_fill_manual(values = c('#4E9E23', '#29BCE3', '#D13D21','grey')) +
  scale_x_continuous(breaks = c(0, 25, 50, 75, 100), limits = c(0, 100)) +
  labs(y = '', x = 'Relative abundance [%]') +
  theme(axis.ticks.y = element_blank(),  axis.text.y = element_blank())
domain

ggarrange(conitgs_domain + labs(tag = 'A. DIAMOND + LCA'), 
          domain + labs(tag = 'B. MetaPhlAn'), ncol = 1, nrow = 2, align = 'v', 
          common.legend = TRUE, 
          legend = 'bottom')
ggsave('out/DIAMOND_metaphlan.svg', dpi=600)
ggsave('out/DIAMOND_metaphlan.png', dpi=600)

# 
contigs %>%
  group_by(Kingdom) %>%
  summarise(no_contigs = n_distinct(contigID)) %>%
  ggplot(aes(x = no_contigs, y = reorder(Kingdom, no_contigs))) +
  geom_col() +
  labs(x = '# contigs', y = '')
ggsave('out/contigs/c_kingdom_contigs.png')

contigs %>%
  group_by(Kingdom) %>%
  summarise(no_contigs = n_distinct(contigID)) %>%
  ggplot(aes(x = no_contigs, y = reorder(Kingdom, no_contigs))) +
  geom_col() +
  labs(x = '# contigs', y = '', subtitle = 'Disparity was less than 1%') +
  scale_x_log10()
ggsave('out/contigs/c_kingdom_disparity_log.png')

# Turn into long file! 
contigs_long <- contigs %>%
  pivot_longer(names_to = 'Group', values_to = 'TPM', cols = starts_with('TPM')) %>%
  mutate(Group = str_remove(Group, 'TPM ')) %>%
  left_join(metadata, by = 'Group') 

# Number of different organisms through time 
c_domain <- contigs_long %>%
  select(contigID, Kingdom, Clade, Phylum, Class, Order, Family, Genus, Species, person, biota, time_point, TPM) %>% 
  group_by(person, biota, time_point, Kingdom) %>%
  reframe(sum = sum(TPM)) 

kingdom_contigs <- c_domain %>% 
  ggplot(aes(x = biota, y = sum, fill = Kingdom)) +
  geom_col() +
  labs(x = '', y = 'TPM', color = 'Kingdom')

# Metaphlan 








c_domain %>%
  ggplot(aes(x = time_point, y = sum, fill = Kingdom)) +
  geom_col() +
  facet_grid(biota~person, scales = 'free_y') +
  labs(x = 'Time point', y = 'TPM', color = 'Individual')
ggsave('out/contigs/c_TPM_kingdom_daily.svg')

##
# Archaea
##
# Filtering big file down for only Archaea
# contigs_archea <- contigs %>%
#   filter(Kingdom == 'Archaea') 
# saveRDS(contigs_archea, 'data/intermediate/conitgs_archaea.RDS')
contigs_archea <- readRDS("~/proj")

# long format 
contigs_archea <- contigs_archea %>%
  pivot_longer(names_to = 'Group', values_to = 'TPM', cols = starts_with('TPM')) %>%
  mutate(Group = str_remove(Group, 'TPM ')) %>%
  left_join(metadata, by = 'Group') 

# Species / Genus  level 
contigs_archea %>%
  filter(biota == 'bulk microbiota') %>%
  ggnested(aes(x = time_point, y = TPM, main_group = Genus, sub_group= Species), 
           legend_title = 'Highest determined taxonomy') +
  geom_col() +
  facet_wrap(~person, nrow = 3, scales = 'free_y') +
  labs(x = 'Time point', subtitle = 'Bulk microbiota samples')
ggsave('plots/contigs/c_archaea_class_species.png')

contigs_archea %>%
  filter(biota == 'ethanol treated sample') %>%
  ggnested(aes(x = time_point, y = TPM, main_group = Genus, sub_group= Species), 
           legend_title = 'Highest determined taxonomy') +
  geom_col() +
  facet_wrap(~person, nrow = 3, scales = 'free_y') +
  labs(x = 'Time point', subtitle = 'Ethanol + EMA treated samples')
ggsave('plots/contigs/c_archaea_class_species_etoh.png')

contigs_archea %>%
  filter(Species == 'Methanobrevibacter smithii') %>%
  group_by(biota, person, time_point) %>%
  summarise(sum = (sum(TPM)/1000000)*100, .groups='drop') %>%
  ggplot(aes(x = time_point, y = sum, color = person)) +
  geom_point(aes(shape=biota), size = 3) +
  geom_line(linewidth = 1) +
  facet_wrap(~biota, scales = 'free') +
  labs(x = 'Time point', y = 'Relative abundance (%)', color = 'Individual', shape = '', subtitle = 'Methanobrevibacter smithii') +
  theme(legend.position = 'bottom')
ggsave('plots/contigs/c_archaea_Msmithii.png')

contigs_archea %>%
  filter(Species == 'Methanosphaera stadtmanae') %>%
  group_by(biota, person, time_point) %>%
  summarise(sum = (sum(TPM)/1000000)*100, .groups='drop') %>%
  ggplot(aes(x = time_point, y = sum, color = person)) +
  geom_point(size = 3) +
  geom_line(linewidth = 1) +
  facet_wrap(~biota, scales = 'free') +
  labs(x = 'Time point', y = 'Relative abundance (%)', color = 'Individual', shape = '', subtitle = 'Methanobrevibacter stadtmanae') +
  theme(legend.position = 'bottom')
ggsave('plots/contigs/c_archaea_Mstadtmanae.png')

# Raw read count table
tab_archaea <- contigs_archea %>% 
  select(contigID, starts_with('Raw read count')) %>% 
  rename_with(~ str_remove(., "^Raw read count "), starts_with("Raw read count")) %>% 
  column_to_rownames('contigID') %>%
  t()

# Alpha diversity 
# Richness

# based on the number of contigs 
alpha <- t(estimateR(tab_archaea)) %>%
  as.data.frame() %>% 
  rownames_to_column('Group') %>%
  left_join(metadata, by = 'Group')

alpha %>%
  ggplot(aes(x = time_point, y = S.obs, color = person)) +
  geom_point(size = 3) +
  geom_line(linewidth = 1) +
  #scale_y_log10() +
  facet_grid(~biota, scales = 'free') +
  labs(x = 'Time point', y = '# archaeal contigs', color = 'Individual', subtitle= '# contigs determined Archaea')
ggsave('plots/contigs/c_archaea_observed.png')

# Alpha diversity based on the summed raw reads mapped to different archaeal genera
genera_archaea <- contigs_archea %>% 
  select(Genus, starts_with('Raw read count')) %>%
  rename_with(~ str_remove(., "^Raw read count "), starts_with("Raw read count")) %>% 
  pivot_longer(names_to = 'Group', values_to = 'read_count', cols = starts_with(c('M', 'S'))) %>%
  group_by(Genus, Group) %>%
  summarise(read_count = sum(read_count), .groups= 'drop') %>%
  pivot_wider(names_from = 'Group', values_from = 'read_count') %>%
  column_to_rownames('Genus') %>%
  t() 

alpha_genera_archea <- genera_archaea %>%
  estimateR() %>%
  t() %>%
  as.data.frame() %>% 
  rownames_to_column('Group') %>%
  left_join(metadata, by = 'Group')

alpha_genera_archea %>%
  ggplot(aes(x = time_point, y = S.obs, color = person)) +
  geom_point(size = 3) +
  geom_line(linewidth = 1) +
  #scale_y_log10() +
  facet_grid(~biota, scales = 'free') +
  labs(x = 'Time point', y = '# archaeal contigs', color = 'Individual', subtitle= '# contigs determined Archaea')

alpha_genera_archea %>%
  ggplot(aes(x = person, y = S.obs, color = person)) +
  geom_boxplot(fill = 'grey', color = 'grey',  alpha = .9, outlier.shape = 4, outlier.size = 2) +
  geom_jitter(size = 2) +
  facet_grid(~biota, scales = 'free') 

# Beta diversity based on genera 
min(rowSums(genera_archaea))

beta_genera_archea <- vegan::avgdist(genera_archaea, sample = 10, iterations = 100, method = 'bray')
nmds <-  metaMDS(beta_genera_archea)
nmds_positions <- as.data.frame(scores(nmds, display='sites')) %>%
  rownames_to_column('Group') %>%
  left_join(metadata %>% select(Group, person, date, biota), by = 'Group')

nmds_positions %>%
  ggplot(aes(x=NMDS1, y=NMDS2, color=person, shape = biota)) +
  geom_jitter(size=4) +
  facet_grid(~biota, scales = 'free') +
  labs(x='', y='', color='Individual')

##
# Eucaryota 
##

# contigs_eukaryota <- contigs %>%
#   filter(Kingdom == 'Eukaryota')
# saveRDS(contigs_eukaryota, 'data/intermediate/contigs_eukaryota.RDS')
contigs_eukaryota <- readRDS("~/projects/longitudinal_shotgun/data/intermediate/contigs_eukaryota.RDS") 

contigs_eukaryota_long <- contigs_eukaryota %>%
  pivot_longer(names_to = 'Group', values_to = 'TPM', cols = starts_with('TPM')) %>%
  mutate(Group = str_remove(Group, 'TPM ')) %>%
  left_join(metadata, by = 'Group')


# Bulk microbiota samples
contigs_eukaryota_filt <- contigs_eukaryota_long %>%
  group_by(person, time_point, biota, Clade, Phylum, Genus) %>%
  reframe(sumTPM = sum(TPM)) %>%
  mutate(Phylum = ifelse(sumTPM < 100, 'Less than 0.01%', Phylum), 
         Genus = ifelse(sumTPM < 100, 'Less than 0.01%', Genus))

contigs_eukaryota_filt %>% 
  filter(biota == 'bulk microbiota') %>%
  ggnested(aes(x = time_point, y = sumTPM, main_group = Clade, sub_group= Phylum), legend_title = '') +
  geom_col() +
  facet_wrap(~person, nrow = 3, scales = 'free_y') +
  labs(x = 'Time point', subtitle = 'Bulk microbiota samples')
ggsave('plots/contigs/c_eucaryota_phylum_bulkmicrobiota_TPMless10.png')

contigs_eukaryota_filt %>%
  filter(biota == 'bulk microbiota') %>%
  ggplot(aes(x = time_point, y = sumTPM, fill = Genus)) +
  geom_col() +
  facet_wrap(~person, nrow = 3, scales = 'free_y') +
  labs(x = 'Time point', subtitle = 'Bulk microbiota samples')
ggsave('plots/contigs/c_eucaryota_genus_bulkmicrobiota_TPMless10.png')


# Ethanol treated samples 
contigs_eukaryota_filt %>%
  filter(biota == 'ethanol treated sample') %>%
  ggplot(aes(x = time_point, y = sumTPM, fill= Phylum)) +
  geom_col() +
  facet_wrap(~person, nrow = 3, scales = 'free_y') +
  labs(x = 'Time point', subtitle = 'Ethanol treated samples')
ggsave('plots/contigs/c_eucaryota_phylum_etoh.png')

contigs_eukaryota_filt  %>%
  filter(biota == 'ethanol treated sample') %>%
  ggplot(aes(x = time_point, y = sumTPM, fill = Genus)) +
  geom_col() +
  facet_wrap(~person, nrow = 3, scales = 'free_y') +
  labs(x = 'Time point', subtitle = 'Ethanol treated samples') 
ggsave('plots/contigs/c_eukaryota_genus_etoh.png')



# Fungi 
contigs_fungi <- contigs_eukaryota_long %>%
  filter(Phylum == 'Fungi') %>%
  group_by(person, time_point, biota, Clade, Phylum, Genus, Species) %>%
  reframe(sumTPM = sum(TPM)) %>%
  mutate(Phylum = ifelse(sumTPM < 100, 'Less than 0.01%', Phylum), 
         Genus = ifelse(sumTPM < 100, 'Less than 0.01%', Genus), 
         Species = ifelse(sumTPM < 100, 'Less than 0.01%', Species))

contigs_fungi %>%
  ggplot(aes(x = time_point, y = (sumTPM/1000000)*100, color = person)) +
  geom_point(size = 3) +
  facet_wrap(~biota, scales = 'free_y') +
  labs(x = 'Time point', y = 'Relative abundance (%)', color = 'Individual', subtitle = 'Relative abundance of fungi')
ggsave('plots/contigs/c_fungi_person_time.png')

# Diversity of fungi in the dataset 

contigs_fungi <- contigs_eukaryota_long %>%
  filter(Phylum == 'Fungi')

tab_fungi <- contigs_eukaryota_long %>%
  filter(Phylum == 'Fungi') %>%
  group_by(Group, Species) %>%
  summarise(TPM = round(sum(TPM), digits = 0), .groups = 'drop') %>%
  pivot_wider(names_from = 'Species', values_from = 'TPM') %>%
  column_to_rownames('Group')

estimateR(tab_fungi) %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column('Group') %>%
  left_join(metadata, by = 'Group') %>%
  ggplot(aes(x = time_point, y = S.obs, color = person)) +
  geom_point(size = 2) +
  geom_line(linewidth = 1) +
  scale_y_continuous(breaks = c(0,2,4,6,8,10)) +
  facet_wrap(~biota) +
  labs(y = '# observed species of fungi', x = 'Time point', color = 'Individual')
ggsave('plots/contigs/c_fungi_species_observed.png')

##
# Virueses 
##
# contigs_viruses <- contigs %>%
#   filter(Kingdom == 'Viruses')
# saveRDS(contigs_viruses, 'data/intermediate/contigs_viruses.RDS')
contigs_viruses <- readRDS("~/projects/longitudinal_shotgun/data/intermediate/contigs_viruses.RDS")

contigs_viruses %>%
  ggplot(aes(x = time_point, y = TPM, fill = Clade)) +
  geom_col() +
  facet_wrap(~person, nrow = 3, scales = 'free') +
  labs(x = 'Time point')
ggsave('plots/contigs/c_viruses_clade.png')

# All viruses were only determined up to Bacteriophage sp. = so I will run geNomad for better classification
# of viral sequences and do the analysis after that

##
# Bacteria 
##

# contigs_bacteria <- contigs %>%
#   filter(Kingdom == 'Bacteria')
# saveRDS(contigs_bacteria, 'data/intermediate/contigs_bacteria.RDS')
contigs_bacteria <- readRDS("~/projects/longitudinal_shotgun/data/intermediate/contigs_bacteria.RDS")

raw_count_bact <- select(contigs_bacteria, Phylum, starts_with('Raw read')) %>%
  pivot_longer(-Phylum) %>%
  group_by(Phylum, name) %>%
  summarise(value = sum(value)) %>%
  mutate(name = str_remove(name, 'Raw read count ')) %>%
  left_join(metadata, by = join_by('name' == 'Group'))

no_reads_sample <- raw_count_bact %>% group_by(name) %>%
  summarise(all_reads = sum(value), .groups = 'drop')

rel_bact <- raw_count_bact %>%
  left_join(no_reads_sample, by = 'name') %>%
  mutate(rel_abund = (value/all_reads)*100)

# Bulk microbiota & ethanol treated samples
rel_bact %>%
  filter(!is.na(Phylum)) %>%
  group_by(biota, person, time_point, Phylum) %>%
  summarise(rel_abund = sum(rel_abund), .groups = 'drop') %>%
  mutate(phylum = ifelse(rel_abund < 1, 'Less than 2 % reads', Phylum)) %>%
  ggplot(aes(x = time_point, y = rel_abund, fill = phylum)) +
  geom_col() +
  facet_grid(biota ~ person)

