# Comparison number of OTUs / Species per genus /per sample etc. 

library(readr)
library(tidyr)
library(dplyr)
library(tibble)
library(vegan)
library(lubridate)
library(ggpubr)
library(purrr)
library(stringr)
library(readxl)

set.seed(96)
theme_set(theme_bw(base_size = 14))

# OTU data 

otutab <- readRDS('data/longitudinal_amplicons/otutab_ethanol_bulk.RDS')
taxtab <- readRDS('data/longitudinal_amplicons/taxtab.RDS')


otu_long <- pivot_longer(as.data.frame(otutab) %>%  rownames_to_column('Group'), cols = starts_with('Otu')) %>%  
  left_join(taxtab, by = 'name')

# number of OTUs per sample
otus <- filter(otu_long, value > 0) %>%  
  group_by(Group) %>% 
  reframe(otus = n_distinct(name))
  
# Number of Genuses per sample 
genus_otu <- filter(otu_long, value > 0) %>%  
  group_by(Group) %>% 
  reframe(genus_otus = n_distinct(Genus))


# metaG data 
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


bacteria <- filter(abund, Domain == 'Bacteria', !is.na(Phylum), !is.na(Class), 
                   !is.na(Order), !is.na(Family), !is.na(Genus), !is.na(Species), !is.na(SGB)) %>% 
  pivot_longer(-c(Domain, Phylum, Class, Order, Family, Genus, Species, SGB)) %>% 
  mutate(PA = ifelse(value > 0, 1, 0), 
         Group = name) %>% 
  select(-name) %>%  
  filter(value > 0)
  

species <- bacteria %>%  
  group_by(Group) %>% 
  reframe(species = n_distinct(Species))

genus_meta <-   bacteria %>%  
  group_by(Group) %>% 
  reframe(genus_meta = n_distinct(Genus))

# Plot together
all <- full_join(otus, species,by = 'Group' ) %>%  
  full_join( genus_meta, by = 'Group') %>% 
  full_join(genus_otu, by = 'Group') %>% 
  left_join(metadata, by = 'Group')


all %>% ggplot(aes(x = otus, y = species, color = biota)) +
  geom_point() +
  geom_smooth(method = 'lm') +
  facet_wrap(~biota, scales = 'free') +
  scale_color_manual(values = c('#3CB371', '#f0a336')) +
  stat_cor(method = "pearson", label.x = 250, label.y = 320) +
  labs(x = '# OTUs/sample ', y = '# species/sample', color = '') +
  theme(legend.position = 'bottom')
ggsave('out/compare_otu_species.png') 

  
all %>% ggplot(aes(x = otus, y = species, color = biota)) +
  geom_point() +
  geom_smooth(method = 'lm') +
  facet_wrap(~biota, scales = 'free') +
  scale_color_manual(values = c('#3CB371', '#f0a336')) +
  stat_cor(method = "pearson", label.x = 250, label.y = 320) +
  labs(x = '# OTUs/sample ', y = '# genera/sample',  color = '') +
  theme(legend.position = 'bottom')
ggsave('out/compare_otu_genera.png')


all %>% ggplot(aes(x = genus_otus, y = genus_meta, color = biota)) +
  geom_point() +
  geom_smooth(method = 'lm') +
  scale_color_manual(values = c('#3CB371', '#f0a336')) +
  facet_wrap(~biota, scales = 'free') +
  stat_cor(method = "pearson") +
  labs(x = '# genera/sample in shotgun metagenomic data', y = '# genera/sample in amplicon data',  color = '') +
  theme(legend.position = 'bottom')
ggsave('out/compare_genera.png')
