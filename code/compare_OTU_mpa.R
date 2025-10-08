# Matching between shotgun and OTUs at the level of genus for bulk microbiota 

library(dplyr)
library(readr)
library(stringr)
library(ggplot2)

theme_set(theme_bw(base_size = 12))

metadata_otu <- readRDS('~/projects/longitudinal_amplicons/data/r_data/metadata.RDS')
metadata_metag <- read.table('~/projects/longitudinal_shotgun/data/metadata.csv', sep = ';', header = T)

otutab <- rownames_to_column(as.data.frame(readRDS('~/projects/longitudinal_amplicons/data/r_data/otutabEM.RDS')), 'Group') %>% 
  pivot_longer(-Group) %>% 
  left_join(readRDS('~/projects/longitudinal_amplicons/data/r_data/taxtab.RDS'), by = 'name') %>% 
  group_by(Group, Domain, Phylum, Class, Order, Family, Genus) %>%
  reframe(value = sum(value)) %>%  
  group_by(Group) %>% 
  mutate(value = value/sum(value)) %>% 
  ungroup()
  

mpatab <- read_tsv('~/projects/longitudinal_shotgun/data/metaphlan_abundance_table.txt', comment = '#') %>%
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
  filter(Domain == 'Bacteria', !is.na(Phylum), !is.na(Class), 
         !is.na(Order), !is.na(Family), !is.na(Genus), !is.na(Species), !is.na(SGB)) %>% 
  pivot_longer(-c(Domain, Phylum, Class, Order, Family, Genus, Species, SGB), names_to = 'Group', values_to = 'value') %>% 
  group_by(Group, Domain, Phylum, Class, Order, Family, Genus) %>%
  reframe(value = sum(value)) %>%  
  filter(Group != 'MC013') %>% 
  left_join(metadata_metag, by = 'Group')

genus_compare <- full_join(otutab, mpatab, by = c('Group', 'Genus')) %>% 
  filter(value.x != 0 | !is.na(value.y)) %>% 
  filter(!is.na(value.x))

genus_compare %>% 
  mutate(diff = abs(value.x - value.y)) %>% 
  ggplot(aes(x = diff, y = Genus)) + 
  geom_point()
