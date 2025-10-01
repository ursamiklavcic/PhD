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

col <- c('#3CB371', '#f0a336')
col2 <- c('#A7E2C1', '#F7CD92')
colm <- '#3CB371'
cole <- '#f0a336'

# Import results from code/sporulation_ability/03_sporeability.sh

blast_results <- read.table('data/sporulation_ability/mpa_blast_results.tsv', sep = '\t', header = TRUE) %>% 
  mutate(genome_id = substr(genome_id, 1, 15)) %>% 
  rename('locus_tag' = 'gene_name')

# Sporulation genes info (weight etc)
gene_info <- read.table('data/sporulation_ability/gene_info.csv', sep = ';', header = TRUE) %>% 
  select(locus_tag, gene_name, weight)

# Import results from 02_find_genomes_from_species.sh
genome_info <- read.table('data/sporulation_ability/mpa_accession_manual.tsv', sep = '\t', header = TRUE)

metadata <- read.table('../longitudinal_shotgun/data/metadata.csv', sep = ';', header = TRUE) %>%  
  mutate(date = dmy(date))


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
sporulation_ability <- spore_ability %>% 
  left_join(genome_info, by = c('genome_id')) 

write.table(sporulation_ability, file='data/sporulation_ability/sporulation_ability2021.tsv',
            quote=FALSE, row.names = FALSE, sep='\t')

### 
# 
abund2 <- read_tsv('~/projects/longitudinal_shotgun/data/metaphlan_abundance_table.txt', comment = '#') %>%
  rename_with(~ str_remove(., '^profiled_'), starts_with('profiled_')) %>%
  filter(grepl('s__', clade_name), !grepl('t__', clade_name)) %>% 
  left_join(select(sporulation_ability, n_genes, PA, sporulation_ability, clade_name), by = 'clade_name') %>% 
  pivot_longer(-c(clade_name, PA, n_genes, sporulation_ability)) %>% 
  #mutate(clade_name = str_remove_all(clade_name, '[a-zA-Z]__')) %>%
  separate(clade_name, into=c('Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species'),
           sep="\\|") %>% 
  mutate(Phylum = ifelse(Phylum == 'p__Firmicutes', 'p__Bacillota', Phylum), 
         Domain = str_remove_all(Domain, 'k__'), 
         Phylum = str_remove_all(Phylum, 'p__'), 
         Class = str_remove_all(Class, 'c__'), 
         Order = str_remove_all(Order, 'o__'), 
         Family = str_remove_all(Family, 'f__'), 
         Genus = str_remove_all(Genus, 'g__'), 
         Species = str_remove_all(Species, 's__')) %>% 
  filter(name != 'MC013') %>% 
  left_join(metadata, by = join_by('name' == 'Group'))

# How many species are undetermined for sporulation ability ? 
length(unique(abund2$Species)) # 1033
length(unique(filter(abund2, is.na(sporulation_ability))$Species)) # 506

abund2 %>% filter(is.na(sporulation_ability), biota == 'bulk microbiota') %>% 
  ggplot(aes(x = as.factor(day), y = value, fill = Phylum)) +
  geom_col() +
  facet_wrap(~person, scales = 'free') +
  labs(x = 'Day', y = 'Relative abundance [%]', subtitle = 'Bulk microbiota')
ggsave('out/sporulation/SGBs_without_sporulation_determination.png')

# SGB/Species that do not have sprulation ability determined represent from sa little
# as 10 up to 35% of all relative abudnance is some samples in bulk microbiota 

abund2 %>% filter(is.na(sporulation_ability), biota == 'ethanol treated sample') %>% 
  ggplot(aes(x = as.factor(day), y = value, fill = Phylum)) +
  geom_col() +
  facet_wrap(~person, scales = 'free') +
  labs(x = 'Day', y = 'Relative abundance [%]', subtitle = 'Ethanol treated samples')
ggsave('out/sporulation/SGBs_without_sporulation_determination_etoh.png')
# And somewhat the same in the ethanol treated samples 

# relative abundance of SGBs with undetermined sporulation ability 
abund2 %>% filter(is.na(sporulation_ability)) %>% 
  ggplot(aes(x = value, y = Phylum, fill = biota)) +
  geom_boxplot() +
  scale_x_log10() +
  scale_fill_manual(values = c('#A7E2C1', '#F7CD92')) +
  labs(x = 'Relative abundance [log10]', y ='', fill = 'Sample type')
ggsave('out/sporulation/relative_abundance_NA.png')

# Prevalence of undefined species 
# plot % species on y, x 0 prevalence % 
preval <- abund2 %>% filter(biota == 'bulk microbiota') %>% 
  group_by(sporulation_ability, person, Phylum, Species) %>% 
  reframe(all_timepoints = n(), 
          timepoints_present = sum(value > 0),
          timepoints_missing = sum(value = 0),
          prevalence = (timepoints_present / all_timepoints) * 100) %>%
  # Calculate number of OTUs per person x treatment x Phylum
  filter(prevalence > 0) %>% 
  group_by(person) %>%
  mutate(all_species_per_person = n_distinct(Species)) %>%
  ungroup() %>%
  group_by(person, sporulation_ability, Phylum, prevalence) %>%
  reframe(n_species_person_spore_phylum = n_distinct(Species), 
          per_species = (n_species_person_spore_phylum / all_species_per_person) * 100) %>% 
  unique() 


preval %>%
  filter(is.na(sporulation_ability)) %>%  
  ggplot(aes(x = prevalence, y = per_species, color = person)) +
  geom_line(linewidth = 2) +
  geom_point(size = 3) +
  # geom_line(data = prevalence %>% filter(sporulation_ability == 'Spore-former'),
  #           aes(group = interaction(person, Phylum)),
  #           color = '#DB674F', linewidth = 0.9, alpha = 0.3) +
  # geom_line(data = prevalence %>% filter(sporulation_ability == 'Non-spore-former'),
  #           aes(group = interaction(person, Phylum)),
  #           color = '#4FC4DB', linewidth = 0.9, alpha = 0.3) +
  # geom_line(data = prevalence %>% filter(is.na(sporulation_ability)),
  #           aes(group = interaction(person, Phylum)),
  #           color = '#7EDB4F', linewidth = 0.9, alpha = 0.3) +
  facet_wrap(~Phylum, nrow = 4, scales = 'free_y') +  
  #scale_color_manual(values = c('#4FC4DB','#DB674F','#7EDB4F')) +
  labs(x = 'Within-individual prevalence\n(% of timepoints present)',
       y = '% of all species within a person', subtitle = 'Bulk microbiota') 
ggsave('out/sporulation/prevalence_NA.png')

# General 
# abund2 %>% filter(biota == 'bulk microbiota') %>% 
abund2 %>% filter(biota == 'ethanol treated sample') %>% 
  group_by(sporulation_ability, person, Phylum, Species) %>% 
  reframe(timepoints_present = sum(value > 0)) %>%
  # Calculate number of OTUs per person x treatment x Phylum
  group_by(person, sporulation_ability, Phylum, timepoints_present) %>%
  reframe(n_species_person_spore_phylum = n_distinct(Species)) %>% 
  unique() %>% 
  filter(is.na(sporulation_ability), timepoints_present > 0) %>% 
  mutate(timepoints_present = ifelse(timepoints_present > 12, 12, timepoints_present), 
         Phylum = ifelse(n_species_person_spore_phylum < 2, '< 2 species in phylum', Phylum)) %>% 
  ggplot(aes(x = as.factor(timepoints_present), y = n_species_person_spore_phylum, fill = Phylum)) +
  geom_col() +
  labs(x = '# timepoints present', y = '# species', caption = 'Number of species present at x timepoints per person, colored by phylum') +
  theme(legend.position = 'bottom')
ggsave('out/sporulation/prevalence_barplot_NA_bulk.png')
ggsave('out/sporulation/prevalence_barplot_NA_etoh.png')


# Exploring the data further! 
abund2 %>% filter(Phylum != c('Bacillota', 'Bacteria_unclassified'), !grepl('SGB', Species), Domain == 'Bacteria') %>% 
  ggplot(aes(x = value, y = Genus, fill = sporulation_ability)) +
  geom_boxplot() +
  scale_x_log10() +
  facet_grid(~sporulation_ability)
ggsave('out/sporulation/relative_abundance_Genus_NA.png', width = 21, height = 29, unit = 'cm')

genus_sporulation <- abund2 %>%  filter(Domain == 'Bacteria') %>% 
  group_by(Family, sporulation_ability) %>% 
  reframe(n = n_distinct(Species)) %>% 
  pivot_wider(names_from = sporulation_ability, values_from = n)

genus_sporulation %>% filter(!is.na(`Spore-former`), !is.na(`Non-spore-former`))

genus_species_both <- abund2 %>% 
  filter(Genus %in% c('Clostridiaceae_unclassified', 'Enterocloster', 
                      'Erysipelatoclostridium', 'Eubacterium', 
                      'Lachnospira', 'Lachnospiraceae_unclassified', 
                      'Ruminococcus')) %>% distinct()

na <- filter(genus_sporulation, is.na(`Spore-former`), is.na(`Non-spore-former`))

# 
# There are some genuses without classification as SF or NSF, but googling them determined whether they were classified as SF in vitro! 

# Also another desision - if most of the species in a genus is SF and there are some SGBs they are also SF! 
# and if most species in genus is NSF, than SGBs are also NSF 
# There are some that have both NSF and SF determination. But are defined as high sporulation signature ( more than 33 sporulation genes), 
# so we will count they are also SF like in the Browne et al 2021 paper! 

# library(diptest)
# sporulation_ability
# 
# family_summary <- abund2 %>%
#   filter(!grepl('FGB', Family) & Domain == 'Bacteria') %>% 
#   group_by(Family) %>% 
#   reframe(mean_score = mean(n_genes, na.rm = TRUE), 
#           sd_score = sd(n_genes, na.rm = TRUE), 
#           modality = ifelse(n_distinct(Species) >= 5 & length(unique(Species)) > 1, dip.test(n_genes)$p.value, 0), 
#           modality = ifelse(modality < 0.5, 'unimodal', 'bimodal'), 
#           PA_family = n_distinct(PA)) %>% 
#   mutate(spore_ability = ifelse(mean_score > 33 & modality == 'unimodal' & PA_family >=1, 'Spore-former', 'Non-spore-former'))

# I have 527 spcies, for which I have determined sporulation ability! I will work with those! 

abund3 <- filter(abund2, !is.na(sporulation_ability), Domain == 'Bacteria', biota == 'bulk microbiota')

# relative abudnance of spore/ non-spore forming 
abund3 %>% 
  ggplot(aes(x = value, y = Phylum, fill = sporulation_ability)) +
  geom_boxplot() +
  scale_x_log10() +
  scale_fill_manual(values = col2) +
  labs(x = 'Relative abundance [log10]', y = '', fill = '') +
  theme(legend.position = 'bottom') +
  stat_compare_means()
ggsave('out/sporulation/SNS_rel_abund.png')

# Is distribution of relative abundances for Bacillota the same for non-spore and spore-forming? 
bacillota <- filter(abund3, Phylum == 'Bacillota')
  
wilcox.test(filter(bacillota, sporulation_ability == 'Spore-former')$value, filter(bacillota, sporulation_ability == 'Non-spore-former')$value, alternative = 'greater')

# Wilcoxon rank sum test with continuity correction
# 
# data:  filter(rel_bacillota, sporulation_ability == "Spore-former")$value and filter(rel_bacillota, sporulation_ability == "Non-spore-former")$value
# W = 257906560, p-value < 2.2e-16
# alternative hypothesis: true location shift is greater than 0

# Number of species 
abund3 %>% 
  mutate(PA = ifelse(value > 0, 1, 0)) %>%  
  filter(PA == 1) %>% 
  group_by(sporulation_ability, Phylum) %>% 
  reframe(sum = n_distinct(Species)) %>% 
  ggplot(aes(x = sum, y = Phylum, fill = sporulation_ability)) +
  geom_col(position = position_dodge(width = 0.9)) +
  geom_text(aes(label = sum),
            position = position_dodge(width = 0.9), hjust = -0.1) +
  scale_fill_manual(values = col2) +
  labs(x = '# Species', y = '', fill = '') +
  theme(legend.position = 'bottom')
ggsave('out/sporulation/SNS_N_species.png')

# density rel abund 
bacillota %>% 
  ggplot(aes(x = value, color = sporulation_ability)) +
  geom_density(linewidth = 2) +
  scale_x_log10() +
  scale_color_manual(values = col2) +
  labs(x = 'Relative abundance [log10]', y = 'Density', color = '')
ggsave('out/sporulation/SNS_density_relabund_bacillota.png')


# Prevalence 
preval3 <- abund3 %>% 
  group_by(sporulation_ability, person, Phylum, Species) %>% 
  reframe(all_timepoints = n(), 
          timepoints_present = sum(value > 0),
          timepoints_missing = sum(value = 0),
          prevalence = (timepoints_present / all_timepoints) * 100) %>%
  # Calculate number of OTUs per person x treatment x Phylum
  filter(prevalence > 0) %>% 
  group_by(person) %>%
  mutate(all_species_per_person = n_distinct(Species)) %>%
  ungroup() %>%
  group_by(person, sporulation_ability, Phylum, prevalence) %>%
  reframe(n_species_person_spore_phylum = n_distinct(Species), 
          per_species = (n_species_person_spore_phylum / all_species_per_person) * 100) %>% 
  unique() 

preval3 %>% 
  filter(Phylum == 'Bacillota') %>% 
  ggplot(aes(x = prevalence, y = per_species, color = person)) +
  geom_line(linewidth = 2) +
  geom_point(size = 2) +
  facet_wrap(~sporulation_ability, nrow = 2, scales = 'free_y') +
  #scale_color_manual(values = col2) +
  labs(x = 'Prevalence [days present within person]', y = expression(paste(italic("Bacillota"), " species [% of species at this prevalence]")), color = '') +
  theme(legend.position = 'bottom')
ggsave('out/sporulation/SNS_prevalence_person_Bacillota.png')

# Statistics for prevalence ? 
wt <- abund3 %>% 
  group_by(sporulation_ability, person, Species) %>% 
  reframe(all_timepoints = n(), 
          timepoints_present = sum(value > 0),
          timepoints_missing = sum(value = 0),
          prevalence = (timepoints_present / all_timepoints) * 100) 

wilcox.test(filter(wt, sporulation_ability == 'Spore-former')$prevalence, filter(wt, sporulation_ability == 'Non-spore-former')$prevalence)

# Wilcoxon rank sum test with continuity correction
# 
# data:  filter(wt, sporulation_ability == "Spore-former")$prevalence and filter(wt, sporulation_ability == "Non-spore-former")$prevalence
# W = 3265983, p-value < 2.2e-16
# alternative hypothesis: true location shift is not equal to 0

# By person 
group_by(wt, person) %>% 
  reframe(wilcox.test(filter(wt, sporulation_ability == 'Spore-former')$prevalence, filter(wt, sporulation_ability == 'Non-spore-former')$prevalence)$p.value)

# 1 A          1.07e-31
# 2 B          1.07e-31
# 3 C          1.07e-31
# 4 D          1.07e-31
# 5 E          1.07e-31
# 6 F          1.07e-31
# 7 G          1.07e-31
# 8 H          1.07e-31
# 9 I          1.07e-31

library(ggpubr)
preval3 %>% 
  filter(Phylum == 'Bacillota') %>% 
  ggplot(aes(x = sporulation_ability, y = prevalence, fill = sporulation_ability)) +
  geom_boxplot() +
  stat_compare_means() +
  scale_fill_manual(values = col2) +
  labs(x = '', y = 'Prevalence [days present within person]') + theme(legend.position = 'none')
ggsave('out/sporulation/SNS_boxplot_prevalence.png')

# Alpha diveristy 
bacillota <- mutate(bacillota, name = paste0(ifelse(sporulation_ability == 'Spore-former', 'S_', 'NS_'),name))

# Alpha diversity 
library(vegan)
n <- filter(bacillota, value > 0) %>% 
  group_by(name) %>% 
  reframe(richness = n_distinct(Species)) 

tab <- bacillota %>% 
  select(Species, name, value) %>%
  pivot_wider(names_from = 'name', values_from = 'value', values_fill = 0) %>% 
  column_to_rownames('Species') %>% 
  t()
shannon = diversity(tab, index = 'shannon')

alpha <- left_join(n, as_tibble(as.list(shannon)) %>% 
                     pivot_longer(names_to = 'name', values_to = 'shannon', cols = starts_with(c('NS', 'S'))), 
                   by = 'name') %>%
  left_join(bacillota, by = 'name') %>% 
  mutate(person2 = person) 

event_data <- read.table('data/extreme_event_data.csv', sep = ',', header = TRUE)

# richness
ggplot(alpha, aes(x=day, y=richness)) +
  geom_rect(data = event_data, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = extremevent_type), inherit.aes = FALSE,
            alpha = 0.6) +
  scale_fill_manual(values = c('white','#d94343', '#d98e43', '#f1f011', '#0c9910', '#3472b7', '#7934b7', '#b73485', '#0f5618')) +
  geom_line(data = dplyr::select(alpha, -person) %>% filter(sporulation_ability == 'Non-spore-former'), 
            aes(group=person2), color= colm, linewidth=0.5, alpha=0.5) +
  geom_line(data=alpha %>% dplyr::select(-person) %>% filter(sporulation_ability == 'Spore-former'), 
            aes(group=person2), color= cole, linewidth=0.5, alpha=0.5) +
  geom_line(data=alpha %>% filter(sporulation_ability == 'Non-spore-former'),
            aes(color=person), color= colm, linewidth=1.2) +
  geom_line(data=alpha %>% filter(sporulation_ability == 'Spore-former'), 
            color=cole, linewidth=1.2) +
  facet_wrap(~person, scales = 'free') +
  labs(x='Day', y= 'Richness', fill = 'Event')
ggsave('out/sporulation/SNS_richness.png')

ggplot(alpha, aes(x=day, y=shannon)) +
  geom_rect(data = event_data, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = extremevent_type), inherit.aes = FALSE,
            alpha = 0.6) +
  scale_fill_manual(values = c('white','#d94343', '#d98e43', '#f1f011', '#0c9910', '#3472b7', '#7934b7', '#b73485', '#0f5618')) +
  geom_line(data = dplyr::select(alpha, -person) %>% filter(sporulation_ability == 'Non-spore-former'), 
            aes(group=person2), color= colm, linewidth=0.5, alpha=0.5) +
  geom_line(data=alpha %>% dplyr::select(-person) %>% filter(sporulation_ability == 'Spore-former'), 
            aes(group=person2), color= cole, linewidth=0.5, alpha=0.5) +
  geom_line(data=alpha %>% filter(sporulation_ability == 'Non-spore-former'),
            aes(color=person), color= colm, linewidth=1.2) +
  geom_line(data=alpha %>% filter(sporulation_ability == 'Spore-former'), 
            color=cole, linewidth=1.2) +
  facet_wrap(~person, scales = 'free') +
  labs(x='Day', y= 'Shannon', fill = 'Event')
ggsave('out/sporulation/SNS_shannon.png')

# Beta diversity & density of clustering

dist <- vegdist(tab, method = 'bray')
mds <- metaMDS(dist)

mds_data <- as.data.frame(mds$points) %>%  
  rownames_to_column('name') %>% 
  mutate(name_old = str_remove_all(name, 'S_'), 
         name_old = str_remove_all(name_old, 'N')) %>% 
  left_join(metadata, by = join_by('name_old' == 'Group')) %>% 
  mutate(spore_ability = substr(name, 1, 2))


ggplot(mds_data, aes(x = MDS1, y = MDS2, color = person, shape = biota)) +
  geom_point()

# 
dist_bray <- as.matrix(dist) %>%
  as_tibble(rownames = 'Group') %>%
  pivot_longer(-Group) %>%
  filter(Group != name) %>%  
  mutate(sample_pairs = paste(Group, name)) %>%
  group_by(sample_pairs) %>%
  reframe(mean_value = mean(value, na.rm = TRUE), 
            median_value = median(value, na.rm = TRUE)) 

# Tidy the Bray data and join with metadata
bray <- dist_bray %>%
  separate(sample_pairs, into = c("Group", "name"), sep = " ") %>%
  mutate(Group_clean = str_remove(Group, "^(NS_|S_)"),
         name_clean = str_remove(name, "^(NS_|S_)"), 
         spore_ability.x = substr(Group, 1, 2), 
         spore_ability.y = substr(name, 1, 2)) %>%
  left_join(metadata %>% select(Group, person, date), by = join_by('Group_clean' == 'Group')) %>%
  left_join(metadata %>% select(Group, person, date), by = join_by('name_clean' == 'Group')) %>%
  mutate(same_person = ifelse(person.x == person.y, 'Within individual', 'Between individuals'),
         same_fraction = ifelse(spore_ability.x == spore_ability.y, 'Yes', 'No')) %>%
  filter(same_fraction == 'Yes')

bray %>%
  mutate(spore_ability.x = ifelse(spore_ability.x == 'S_', 'Spore-forming', 'Non-spore-forming')) %>% 
  ggplot(aes(x=spore_ability.x, y=median_value, fill=spore_ability.y)) +
  geom_boxplot() +
  stat_compare_means() +
  scale_fill_manual(values = col2) +
  facet_wrap(~same_person) +
  labs(x = 'Community', y = 'Median Bray-Curtis dissimilarity') +
  theme(legend.position = 'none')
ggsave('out/sporulation/SNS_bray_curtis_boxplot.png')

# in time 
time_bray <-  as.matrix(dist) %>%
  as_tibble(rownames = 'Group') %>%
  pivot_longer(-Group) %>%
  filter(Group != name) %>%  
  mutate(Group_clean = str_remove(Group, "^(NS_|S_)"),
         name_clean = str_remove(name, "^(NS_|S_)"), 
         spore_ability.x = substr(Group, 1, 2), 
         spore_ability.y = substr(name, 1, 2)) %>% 
  left_join(metadata %>% select(Group, person, date), by = join_by('Group_clean' == 'Group')) %>%
  left_join(metadata %>% select(Group, person, date), by = join_by('name_clean' == 'Group')) %>%
  mutate(same_person = ifelse(person.x == person.y, 'Within individual', 'Between individuals'),
         same_fraction = ifelse(spore_ability.x == spore_ability.y, 'Yes', 'No'),
         spore_ability.x = ifelse(spore_ability.x == 'S_', 'Spore-forming', 'Non-spore-forming')) %>%
  filter(same_fraction == 'Yes') %>%  
  # Filter different individuals
  filter(same_person == 'Within individual') %>%
  # Calculate the difference between sampling times
  mutate(diff=abs(date.x-date.y)) %>%
  # group by difference between days and person
  group_by(spore_ability.x, person.x, diff) %>%
  reframe(median=median(value)) 

time_bray %>%
  ggplot(aes(x=diff, y=median, color=spore_ability.x)) +
  geom_point() +
  geom_smooth(method = 'lm') +
  stat_cor(method = 'pearson', aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), size = 4) +
  labs(x='Days between sampling', y='Median Bray-Curtis distance', color='') +
  theme(legend.position = 'bottom') +
  guides(fill = guide_legend(ncol = 2)) +
  scale_color_manual(values = col2) 
ggsave('out/sporulation/SNS_time_braycurtis.png')

