library(stringr)
library(ggplot2)
library(readr)
library(dplyr)
library(tibble)
library(lubridate)
library(tidyr)
library(ggnewscale)
library(vegan)

theme_set(theme_bw())
col <- c('#3CB371', '#f0a336')
col2 <- c('#A7E2C1', '#F7CD92')
colm <- '#3CB371'
cole <- '#f0a336'
# # Windows 
# metadata <- read_csv2("G:/projekti/longitudinal_shotgun/data/metadata.csv")
# 
# abund <- read_tsv("G:/projekti/longitudinal_shotgun/data/metaphlan_abundance_table.txt", comment = '#') %>%
#   rename_with(~ str_remove(., '^profiled_'), starts_with('profiled_')) %>%
#   mutate(clade_name = str_remove_all(clade_name, '[a-zA-Z]__')) %>%
#   separate(clade_name, into=c('Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species'),
#            sep="\\|")
#   # mutate(Clade = ifelse(is.na(Clade), Kingdom, Clade), 
#   #        Phylum = ifelse(is.na(Phylum), Clade, Phylum), 
#   #        Class = ifelse(is.na(Class), Phylum, Class), 
#   #        Order = ifelse(is.na(Order), Class, Order), 
#   #        Family = ifelse(is.na(Family), Order, Family), 
#   #        Genus = ifelse(is.na(Genus), Family, Genus), 
#   #        Species = ifelse(is.na(Species), Genus, Species)) %>%
#   

# HPC 
metadata <- read.table('~/projects/longitudinal_shotgun/data/metadata.csv', header= TRUE, sep = ';') %>%
  mutate(date = dmy(date))

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

length(unique(filter(abund, Domain == 'Bacteria')$Species))

abund_gtdb <- read.table('~/projects/longitudinal_shotgun/data/gtdb_merged.txt', sep = '\t', header = TRUE) %>% 
  pivot_longer(-clade_name) %>% 
  filter(grepl('s__', clade_name)) %>% 
  separate(clade_name, into=c('Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species'),
           sep=";")

length(unique(filter(abund_gtdb, Domain == 'd__Bacteria')$Species))

# At kingdom level
domain <- filter(abund, is.na(Phylum)) %>% 
  select(-c(Phylum, Class, Order, Family, Genus, Species, SGB)) %>% 
  pivot_longer(-Domain) %>%  
  left_join(metadata, by = join_by('name' == 'Group'))

domain %>% 
  group_by(biota) %>% 
  mutate(rel_abund = value /sum(value) * 100) %>% 
  ggplot(aes(x = biota, y = rel_abund, fill = Domain)) +
  geom_col() +
  labs(x = '', y = 'Relative abundance [%]')
ggsave('longitudinal_shotgun/plots/mpa_rel_abund_kingdom.png')    

# Phylum level 
phylum <- filter(abund, is.na(Class), !is.na(Phylum), Domain == 'Bacteria') %>% 
  select(-c(Domain, Class, Order, Family, Genus, Species, SGB)) %>% 
  pivot_longer(-Phylum) %>%  
  left_join(metadata, by = join_by('name' == 'Group')) 

phylum %>%  
  group_by(biota) %>% 
  mutate(rel_abund = value /sum(value) * 100) %>% 
  ggplot(aes(x = biota, y = rel_abund, fill = Phylum)) +
  geom_col() +
  labs(x = '', y = 'Relative abundance [%]')
ggsave('longitudinal_shotgun/plots/mpa_rel_abund_phylum.png') 


# how many species did we recover in each phylum 
bacteria <- filter(abund, Domain == 'Bacteria', !is.na(Phylum), !is.na(Class), 
                   !is.na(Order), !is.na(Family), !is.na(Genus), !is.na(Species), !is.na(SGB)) %>% 
  pivot_longer(-c(Domain, Phylum, Class, Order, Family, Genus, Species, SGB)) %>% 
  left_join(metadata, by = join_by('name' == 'Group')) %>% 
  mutate(PA = ifelse(value > 0, 1, 0))

bacteria %>% 
  ggplot(aes(x = value, y = Phylum, fill = biota)) +
  geom_boxplot() + 
  scale_x_log10() +
  scale_fill_manual(values = col2) +
  labs(x = 'Relative abundance [log10]', y = '', fill = 'Sample type') +
  theme(legend.position = 'bottom')
ggsave('longitudinal_shotgun/plots/rel_abund_phylum_boxplot.png')

bacteria %>%
  filter(PA == 1) %>% 
  group_by(biota, Phylum) %>% 
  reframe(sum = n_distinct(Species)) %>% 
  ggplot(aes(x = sum, y = Phylum, fill = biota)) +
  geom_col(position = position_dodge(width = 0.9)) +
  geom_text(aes(label = sum),
    position = position_dodge(width = 0.9), hjust = -0.1) +
  scale_fill_manual(values = col2) +
  labs(x = '# Species', y = '', fill = 'Sample\ntype') +
  theme(legend.position = 'bottom')
ggsave('longitudinal_shotgun/plots/n_species_phylum.png')

# Archea
archaea <- filter(abund, Domain == 'Archaea', !is.na(Species), !is.na(SGB)) %>% 
  pivot_longer(-c(Domain, Phylum, Class, Order, Family, Genus, Species, SGB)) %>%  
  left_join(metadata, by = join_by('name' == 'Group')) 

archaea %>% 
  filter(value > 0) %>%  
  ggplot(aes(x = as.factor(day), y = value, fill = Species)) +
  geom_col(position = 'dodge') +
  facet_wrap(~person, scales = 'free', nrow = 5) +
  scale_fill_manual(values = col2) +
  labs(x = 'Day', y = 'Relative abundance [%]') +
  theme(legend.position = 'bottom') 
ggsave('longitudinal_shotgun/plots/archaea_species.png')

archaea %>% 
  filter(value > 0) %>%  
  mutate(species_SGB = paste(Species, SGB)) %>% 
  ggplot(aes(x = day, y = value, color = species_SGB)) +
  geom_point(size = 2) +
  geom_line(linewidth = 1.5) +
  facet_wrap(~person, scales = 'free', nrow = 5) +
  labs(x = 'Day', y = 'Relative abundance [%]', color = 'Species \n&\nSGB') +
  theme(legend.position = 'bottom') +
  guides(color = guide_legend(nrow = 3, byrow = TRUE))
ggsave('longitudinal_shotgun/plots/archaea_SGB.png')

# Eukaryota
eukaryota <- filter(abund, Domain == 'Eukaryota', !is.na(Species), !is.na(SGB)) %>% 
  pivot_longer(-c(Domain, Phylum, Class, Order, Family, Genus, Species, SGB)) %>%  
  left_join(metadata, by = join_by('name' == 'Group')) 

filter(eukaryota, value > 0) %>%  
  mutate(species_SGB = paste(Species, SGB)) %>% 
  ggplot(aes(x = day, y = value, color = species_SGB)) +
  geom_point(size = 2) +
  geom_line(linewidth = 1.5) +
  facet_wrap(~person, scales = 'free', nrow = 5) +
  labs(x = 'Day', y = 'Relative abundance [%]', color = 'Species \n&\nSGB') +
  theme(legend.position = 'bottom') +
  guides(color = guide_legend(nrow = 2, byrow = TRUE))
ggsave('longitudinal_shotgun/plots/eukaryota_SGB.png')

# Alpha diversity 
n <- filter(bacteria, value > 0) %>% 
  group_by(name) %>% 
  reframe(richness = n_distinct(SGB))

tab <- filter(abund, Domain == 'Bacteria', !is.na(Phylum), !is.na(Class), 
              !is.na(Order), !is.na(Family), !is.na(Genus), !is.na(Species), !is.na(SGB)) %>% 
  select(-c(Domain, Phylum, Class, Order, Family, Genus, Species)) %>% 
  column_to_rownames('SGB') %>% 
  t()
shannon = diversity(tab, index = 'shannon')

alpha <- left_join(n, as_tibble(as.list(shannon)) %>% 
                     pivot_longer(names_to = 'name', values_to = 'shannon', cols = starts_with(c('M', 'S'))), 
                   by = 'name') %>%
  left_join(metadata, by = join_by('name' == 'Group')) %>% 
  mutate(person2 = person) 

# event data
# event_data <- metadata %>%
#   select(person, day, extremevent_type) %>%
#   distinct() %>%
#   filter(!is.na(day)) %>% 
#   mutate(xmin = day - 2, xmax = day +2, ymin = -Inf,ymax = Inf) 
# write_csv(event_data, 'extreme_event_data.csv')
# correct the antibiotics data for person H 

event_data <- read.table('data/extreme_event_data.csv', sep = ',', header = TRUE)

# richness
ggplot(alpha, aes(x=day, y=richness)) +
  geom_rect(data = event_data, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = extremevent_type), inherit.aes = FALSE,
            alpha = 0.6) +
  scale_fill_manual(values = c('white','#d94343', '#d98e43', '#f1f011', '#0c9910', '#3472b7', '#7934b7', '#b73485', '#0f5618')) +
  geom_line(data = dplyr::select(alpha, -person) %>% filter(biota == 'bulk microbiota'), 
            aes(group=person2), color= colm, linewidth=0.5, alpha=0.5) +
  geom_line(data=alpha %>% dplyr::select(-person) %>% filter(biota == 'ethanol treated sample'), 
            aes(group=person2), color= cole, linewidth=0.5, alpha=0.5)+
  geom_line(data=alpha %>% filter(biota == 'bulk microbiota'),
            aes(color=person), color= colm, linewidth=1.2) +
  geom_line(data=alpha %>% filter(biota == 'ethanol treated sample'), 
            color=cole, linewidth=1.2) +
  facet_wrap(~person, scales = 'free') +
  labs(x='Day', y= 'Richness', fill = 'Event')
ggsave('longitudinal_shotgun/plots/mpa_richness.png')

# Shannon 
ggplot(alpha, aes(x=day, y=shannon)) +
  geom_rect(data = event_data, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = extremevent_type), inherit.aes = FALSE,
            alpha = 0.6) +
  scale_fill_manual(values = c('white','#d94343', '#d98e43', '#f1f011', '#0c9910', '#3472b7', '#7934b7', '#b73485', '#0f5618')) +
  geom_line(data = dplyr::select(alpha, -person) %>% filter(biota == 'bulk microbiota'), 
            aes(group=person2), color= colm, linewidth=0.5, alpha=0.5) +
  geom_line(data=alpha %>% dplyr::select(-person) %>% filter(biota == 'ethanol treated sample'), 
            aes(group=person2), color= cole, linewidth=0.5, alpha=0.5)+
  geom_line(data=alpha %>% filter(biota == 'bulk microbiota'),
            aes(color=person), color= colm, linewidth=1.2) +
  geom_line(data=alpha %>% filter(biota == 'ethanol treated sample'), 
            color=cole, linewidth=1.2) +
  facet_wrap(~person, scales = 'free') +
  labs(x='Day', y= 'Shannon', fill = 'Event')
ggsave('longitudinal_shotgun/plots/mpa_shannon.png')

# Composition
bacteria %>% 
  filter(biota == 'ethanol treated sample') %>%
  ggplot(aes(x = factor(day), y = value, fill = Phylum)) +
  geom_col() +
  facet_wrap(~person, scales = 'free_x') +
  labs(x = 'Day', y = 'Relative abundance [%]')
ggsave('longitudinal_shotgun/plots/mpa_rel_abund_bact_etoh.png')

# Sporulation ability 
# spores <- read.table('longitudinal_shotgun/data/sporulation_ability2021.tsv', sep = '\t', header = TRUE)
# 
# # GTDB 
# # joining by whole taxonomical name is inefficient, so I will do it by hand, 
# # to be able to compare metaphlan results with sporulation ability and later PStrain with sporulation ability 
# n <- read.table('~/projects/longitudinal_shotgun/data/gtdb_merged.txt', sep = '\t', header = TRUE) %>% 
#   pivot_longer(-clade_name) %>% 
#   filter(grepl('s__', clade_name)) %>% 
#   left_join(g2, by = join_by('clade_name' == 'tax')) %>% 
#   select(clade_name, genome_id) %>% unique()
# 
# write_tsv(n, 'longitudinal_shotgun/data/metaphlan_gtdb_species_pre.tsv')
# 
# # Import metaphlan_gtdb_species.tsv (which was manually curated) 
# 
# read.table('longitudinal_shotgun/data/metaphlan_gtdb_species.tsv') %>% 
#   left_join(spores, by = join_by('clade_name' == 'tax'), relationship = "many-to-many") %>% 
#   separate(clade_name, into=c('Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species'), sep=";") %>% 
#   mutate(PA = ifelse(value > 0, 1, 0), 
#          name = str_remove_all(name, 'gtdbdd_'), 
#          name = str_remove_all(name, '.txt'),
#          Domain = str_remove_all(Domain, 'd__'), 
#          Phylum = str_remove_all(Phylum, 'p__'), 
#          Class = str_remove_all(Class, 'c__'), 
#          Order = str_remove_all(Order, 'o__'), 
#          Family = str_remove_all(Family, 'f__'), 
#          Genus = str_remove_all(Genus, 'g__'), 
#          Species = str_remove_all(Species, 's__') 
#          # Phylum = str_replace_all(Phylum, 'Firmicutes', 'Bacillota'), 
#          # Phylum = str_replace_all(Phylum, ' Bacillota', 'Bacillota'),
#          # Phylum = str_remove_all(Phylum, '_[a-zA-Z]'),
#          # Phylum = str_remove_all(Phylum, ' '), 
#          # Phylum = str_replace_all(Phylum, 'Actinobacteriota', 'Actinomycetota'),
#          # Phylum = str_replace_all(Phylum, 'Cyanobacteria', 'Cyanobacteriota'), 
#          # Phylum = str_replace_all(Phylum, 'Desulfobacterota_I', 'Desulfobacterota'), 
#          # Phylum = str_replace_all(Phylum, 'Proteobacteria', 'Pseudomonadota'), 
#          # Phylum = if_else(is.na(Phylum), 'unclassified Bacteria', Phylum)
#          ) %>% 
#   filter(Domain == 'Bacteria', !is.na(Species)) %>% 
#   left_join(metadata, by = join_by('name' == 'Group'), relationship = "many-to-many") 
# 
# abund_gtdb %>% 
#   group_by(sporulation_ability) %>% 
#   reframe(n = n_distinct(Species))
# 
