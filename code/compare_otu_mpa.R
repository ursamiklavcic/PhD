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
  left_join(metadata, by = join_by('name' == 'Group')) %>% 
  mutate(PA = ifelse(value > 0, 1, 0))

species <- bacteria %>%
  filter(value > 0) %>% 
  group_by(name) %>% 
  reframe(species = n_distinct(Species))

genus_meta <- bacteria %>% 
  filter(value > 0) %>% 
  group_by(name) %>% 
  reframe(genus_meta = n_distinct(Genus))

# Plot together
all <- full_join(otus, species,by = join_by('Group' == 'name')) %>%  
  full_join(genus_meta, by = join_by('Group' == 'name')) %>% 
  full_join(genus_otu, by = 'Group') %>% 
  left_join(metadata, by = 'Group') %>% 
  filter(!is.na(biota))


all %>% ggplot(aes(x = otus, y = species, color = biota)) +
  geom_point(size = 3) +
  #geom_smooth(method = 'lm') +
  facet_wrap(~biota, scales = 'free') +
  scale_color_manual(values = c('#3CB371', '#f0a336')) +
  geom_abline() +
  labs(x = '# OTUs/sample ', y = '# species/sample', color = '') +
  theme(legend.position = 'none')
ggsave('out/compare_otu_species.png') 

  
all %>% ggplot(aes(x = otus, y = species, color = biota)) +
  geom_point() +
  facet_wrap(~biota, scales = 'free') +
  scale_color_manual(values = c('#3CB371', '#f0a336')) +
  geom_abline() +
  labs(x = '# OTUs/sample ', y = '# genera/sample',  color = '') +
  theme(legend.position = 'none')
ggsave('out/compare_otu_genera.png')


all %>% ggplot(aes(x = genus_otus, y = genus_meta, color = biota)) +
  geom_point() +
  scale_color_manual(values = c('#3CB371', '#f0a336')) +
  facet_wrap(~biota, scales = 'free') +
  geom_abline()+  
  labs(x = '# genera/sample in shotgun metagenomic data', y = '# genera/sample in amplicon data',  color = '') +
  theme(legend.position = 'none')
ggsave('out/compare_genera.png')

# Compare relative abundance 
# col_phylum = c('#1F77B4', '#FF7F0E',  '#2CA02C',  '#D62728', '#9467BD', '#8C564B', '#f4d03f', 
#                )
# otu_rel <- otu_long %>% 
#   left_join(metadata, by = 'Group') %>% 
#   group_by(biota) %>% 
#   mutate(rel_abund = value/sum(value)*100) %>% 
#   ungroup() %>%
#   group_by(biota, Phylum) %>% 
#   reframe(rel_abund_otu = sum(rel_abund)) %>% 
#   mutate(Phylum = case_when(
#     Phylum == 'TM7' ~ 'Saccharibacteria', 
#     Phylum == 'Tenericutes' ~ 'Mycoplasmatota', 
#     Phylum == 'Lentisphaerae' ~ 'Lentisphaerota', 
#     Phylum == 'Synergistetes' ~ 'Synergistota',
#     Phylum == 'Deferribacteres' ~ 'Deferribacteraceae', 
#     Phylum == 'Fusobacteria' ~ 'Fusobacterium',
#     TRUE ~ Phylum))
#  
# unique(otu_rel$Phylum)
# 
# mpa_rel <- bacteria %>% 
#   mutate(value = ifelse(is.na(value), 0, value)) %>% 
#   group_by(biota, Phylum, Species) %>% 
#   reframe(value = value/sum(value, na.rm = TRUE)) %>% 
#   group_by(biota, Phylum) %>% 
#   reframe(rel_abund_mpa = sum(value, na.rm = TRUE)) %>% 
#   mutate(Phylum = case_when(
#     Phylum == 'Actinobacteria' ~ 'Actinomycetota',
#     Phylum == 'Candidatus_Saccharibacteria' ~ 'Saccharibacteria',
#     Phylum == 'Chloroflexi' ~ 'Chloroflexota',
#     Phylum == 'Proteobacteria' ~ 'Pseudomonadota',
#     Phylum == 'Lentisphaerae' ~ 'Lentisphaerota', 
#     Phylum == 'Fusobacteria' ~ 'Fusobacterium',
#     Phylum == 'Verrucomicrobia' ~ 'Verrucomicrobiota', 
#     Phylum == 'Bacteria_unclassified' ~ 'unclassified Bacteria',
#     Phylum == 'Candidatus_Melainabacteria' ~ 'Cyanobacteria', 
#     Phylum == 'Synergistetes' ~ 'Synergistota',
#     Phylum == 'Tenericutes' ~ 'Mycoplasmatota', 
#     TRUE ~ Phylum ))
# 
# unique(mpa_rel$Phylum)
# 
# rel <- full_join(mpa_rel, otu_rel, by = c('biota', 'Phylum')) %>% 
#   pivot_longer(names_to = 'sample', values_to = 'rel_abund', cols = starts_with('rel_abund')) 
# 
# rel %>% 
#   ggplot(aes(x = biota, y = rel_abund, fill = Phylum)) +
#   geom_col() +
#   facet_wrap(~sample)

# Alpha diveristy 

alpha_mpa <- readRDS('data/longitudinal_shotgun/alpha_diveristy.RDS') %>% select(name, richness, shannon)
alpha_otu <- readRDS('data/longitudinal_amplicons/alpha.RDS')

alpha <- full_join(alpha_otu, alpha_mpa, by = join_by('Group' == 'name'))

alpha %>% filter(biota == 'Bulk microbiota') %>% 
  pivot_longer(names_to = 'sample', values_to = 'shannon', cols = starts_with('shannon')) %>% 
  ggplot(aes(x=day, y=shannon)) +
  geom_line(data=alpha_meta %>% dplyr::select(-person) %>% filter(biota == 'Bulk microbiota'), 
            aes(group=person2), color= colm, linewidth=0.5, alpha=0.5) +
  geom_line(data=alpha_meta %>% dplyr::select(-person) %>% filter(biota == 'Ethanol treated sample'), 
            aes(group=person2), color= cole, linewidth=0.5, alpha=0.5)+
  geom_line(data=alpha_meta %>% filter(biota == 'Bulk microbiota'),
            aes(color=person), color= colm, linewidth=1.2) +
  geom_line(data=alpha_meta %>% filter(biota == 'Ethanol treated sample'), 
            color=cole, linewidth=1.2) +
  facet_wrap(~person, scales = 'free') +
  labs(x='Day', y= '#OTUs', color = 'Sample'))
