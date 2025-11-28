library(stringr)
library(ggplot2)
library(readr)
library(dplyr)
library(tibble)
library(lubridate)
library(tidyr)
library(ggnewscale)
library(vegan)
library(ggpubr)

theme_set(theme_bw(base_size=14))
col <- c('#3CB371', '#f0a336')
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

# abund_gtdb <- read.table('~/projects/longitudinal_shotgun/data/gtdb_merged.txt', sep = '\t', header = TRUE) %>% 
#   pivot_longer(-clade_name) %>% 
#   filter(grepl('s__', clade_name)) %>% 
#   separate(clade_name, into=c('Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species'),
#            sep=";")
# 
# length(unique(filter(abund_gtdb, Domain == 'd__Bacteria')$Species))

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
ggsave('out/metaphlan/mpa_rel_abund_kingdom.png')    

# Archea
archaea <- filter(abund, Domain == 'Archaea', !is.na(Species), !is.na(SGB)) %>% 
  pivot_longer(-c(Domain, Phylum, Class, Order, Family, Genus, Species, SGB)) %>%  
  left_join(metadata, by = join_by('name' == 'Group')) 

ag <- archaea %>% 
  filter(value > 0) %>%  
  ggplot(aes(x = person, y = value/n_distinct(name), fill = Genus)) +
  geom_col() +
  scale_fill_manual(values = c('#27CFF5', '#671BA6')) +
  labs(x ='Individual', y = 'Mean relative abundance across samples of a person [%]') +
  theme(legend.position = 'bottom')
ag
ggsave('out/metaphlan/archaea_mean_rel_across.png', dpi = 600)

as <- archaea %>% 
  filter(value > 0) %>%  
  ggplot(aes(x = as.factor(day), y = value, fill = Species)) +
  geom_col(position = 'dodge') +
  scale_fill_manual(values = c('#27CFF5', '#671BA6')) +
  facet_wrap(~person, scales = 'free', nrow = 3) +
  labs(x = 'Day', y = 'Relative abundance [%]') +
  theme(legend.position = 'bottom') 
as
ggsave('out/metaphlan/archaea_species.png', dpi=600)
ggsave('out/metaphlan/archaea_species.svg')

ggarrange(ag + labs(tag = 'A'), as + labs(tag = 'B'), 
          widths = c(0.6, 1))

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
ggsave('out/metaphlan/archaea_SGB.png')

# Eukaryota
eukaryota <- filter(abund, Domain == 'Eukaryota', !is.na(Species), !is.na(SGB)) %>% 
  pivot_longer(-c(Domain, Phylum, Class, Order, Family, Genus, Species, SGB)) %>%  
  left_join(metadata, by = join_by('name' == 'Group')) 

filter(eukaryota, value > 0) %>%  
  mutate(#species_SGB = paste(Species, SGB), 
    species = ifelse(str_detect(Species, 'Blastocystis'), 'Blastocystis sp.', 'Saccharomyces cerevisiae')) %>% 
  ggplot(aes(x = as.factor(day), y = value, fill = species)) +
  geom_col(position = 'dodge') +
  scale_fill_manual(values = c('#2761F5', '#4BC4BF')) +
  facet_wrap(~person, scales = 'free', nrow = 3) +
  labs(x = 'Day', y = 'Relative abundance [%]', fill = 'Species') +
  theme(legend.position = 'bottom') +
  guides(color = guide_legend(nrow = 2, byrow = TRUE))
ggsave('out/metaphlan/eukaryota_SGB.png', dpi = 600)
ggsave('out/metaphlan/eukaryota_SGB.svg')

# Bacteria 
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
ggsave('out/metaphlan/mpa_rel_abund_phylum.png', dpi = 600) 

# how many species did we recover in each phylum 
bacteria <- filter(abund, Domain == 'Bacteria', !is.na(Phylum), !is.na(Class), 
                   !is.na(Order), !is.na(Family), !is.na(Genus), !is.na(Species), !is.na(SGB)) %>% 
  pivot_longer(-c(Domain, Phylum, Class, Order, Family, Genus, Species, SGB)) %>% 
  left_join(metadata, by = join_by('name' == 'Group')) %>% 
  mutate(PA = ifelse(value > 0, 1, 0)) %>% 
  filter(!is.na(biota))

bacteria %>% 
  ggplot(aes(x = value, y = Phylum, fill = biota)) +
  geom_boxplot() + 
  scale_x_log10() +
  scale_fill_manual(values = col) +
  labs(x = 'Relative abundance [log10]', y = '', fill = '') +
  theme(legend.position = 'bottom')
ggsave('out/metaphlan/rel_abund_phylum_boxplot.png')

# For together
rel <- bacteria %>%
  ggplot(aes(x = value, y = Phylum, fill = biota)) +
  geom_boxplot() + 
  scale_x_log10() +
  scale_fill_manual(values = col) +
  labs(x = 'Relative abundance [log10]', y = '', fill = '') +
  theme(legend.position = 'bottom',  axis.text.y=element_blank(), axis.ticks.y=element_blank())
rel

no <- bacteria %>%
  filter(PA == 1) %>% 
  group_by(biota, Phylum) %>% 
  reframe(sum = n_distinct(Species)) %>% 
  ggplot(aes(x = sum, y = Phylum, fill = biota)) +
  geom_col(position = position_dodge(width = 0.9)) +
  geom_text(aes(label = sum),
    position = position_dodge(width = 0.9), hjust = -0.1) +
  scale_fill_manual(values = col) +
  labs(x = '# Species', y = '', fill = '') +
  theme(legend.position = 'bottom')
no
ggsave('out/metaphlan/n_species_phylum.png')

# both 
ggarrange(no + labs(tag = 'A'), rel + labs(tag = 'B'), 
          common.legend = TRUE, legend = 'bottom', 
          widths = c(1, .7))
ggsave('out/metaphlan/mpa_number_species_relabund.png')

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
alpha
saveRDS(alpha, 'data/longitudinal_shotgun/alpha_diveristy.RDS')

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
ggsave('out/metaphlan/mpa_richness.png')

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
ggsave('out/metaphlan/mpa_shannon.png')

# Composition
bacteria %>% 
  filter(biota == 'Bulk microbiota') %>%
  ggplot(aes(x = factor(day), y = value, fill = Phylum)) +
  geom_col() +
  facet_wrap(~person, scales = 'free_x') +
  labs(x = 'Day', y = 'Relative abundance [%]')
ggsave('out/metaphlan/mpa_rel_abund_bact.png')

bacteria %>% 
  filter(biota == 'ethanol treated sample') %>%
  ggplot(aes(x = factor(day), y = value, fill = Phylum)) +
  geom_col() +
  facet_wrap(~person, scales = 'free_x') +
  labs(x = 'Day', y = 'Relative abundance [%]')
ggsave('out/metaphlan/mpa_rel_abund_bact_etoh.png')


# Beta diveristy 
tab <- filter(bacteria, biota == 'bulk microbiota') %>% select(SGB, name, value) %>% 
  pivot_wider(names_from = 'name', values_from = 'value', values_fill = 0) %>% 
  column_to_rownames('SGB') %>%  
  t()

dist <- vegdist(tab, method = 'bray')
saveRDS(dist, 'data/longitudinal_shotgun/dist_mpa.RDS')
nmds <- metaMDS(dist)

nmds_positions <- as.data.frame(scores(nmds, display="sites")) %>%
  rownames_to_column('Group')

dist_meta = nmds_positions %>%
  left_join(metadata, by = 'Group')
saveRDS(dist_meta, 'data/longitudinal_shotgun/nmds_mpa_positions.rds')

dist_meta %>%
  ggplot(aes(x=NMDS1, y=NMDS2, color=person, shape = biota)) +
  geom_point(size = 5) +
  labs(color = 'Individual', shape = '') +
  theme(legend.position = 'bottom')
ggsave('out/metaphlan/mps_nmds_bray_sample_type.png')

# For only microbiota 
nmds_mpa <- dist_meta %>% filter(biota == 'bulk microbiota') %>% 
  ggplot(aes(x=NMDS1, y=NMDS2, color=person)) +
  geom_point(size = 5) +
  labs(color = 'Individual') 

nmds_otu <- readRDS('out/longitudinal_amplicons/nmds_otu.RDS')

ggarrange(nmds_otu +labs(tag = 'A'), nmds_mpa +labs(tag = 'B'), common.legend = T, 
          legend = 'bottom')
ggsave('out/nmds_both.png', dpi = 600)

dist_meta %>%
  mutate(extreme = ifelse(extremevent_type != '', 'extreme event', '')) %>% 
  ggplot(aes(x=NMDS1, y=NMDS2, color=person, shape = extreme)) +
  geom_point(size = 5) +
  labs(color = 'Individual', shape = '') +
  theme(legend.position = 'bottom')
ggsave('out/metaphlan/mpa_nmds_bray_extreme.png')

# Distances between individuals samples and extreme events  etc. 
# Tidy the Bray-Curtis matrix
dist_long <- as.matrix(dist) %>%
  as_tibble(rownames = 'Group') %>%
  pivot_longer(-Group) %>%
  filter(Group != name)

meta <- select(metadata, Group, person, date, biota)

dist_all <- dist_long %>%
  mutate(sample_pairs = paste(Group, name)) %>%
  group_by(sample_pairs) %>%
  reframe(mean_value = mean(value, na.rm = TRUE), 
            median_value = median(value, na.rm = TRUE)) %>%
  separate(sample_pairs, into = c("Group", "name"), sep = " ") %>%
  left_join(meta, by = 'Group') %>%
  left_join(meta, by = join_by('name' == 'Group')) %>%
  mutate(same_person = ifelse(person.x == person.y, 'Intra individual', 'Inter individual'), 
         same_sample = ifelse(biota.x == biota.y, 'Same', 'Different'), 
         date_dist = abs(date.x-date.y))

dist_all %>%  filter(same_sample == 'Same') %>%  
  ggplot() +
  geom_boxplot(mapping = aes(x = same_person, y = median_value, fill = biota.x)) +
  scale_fill_manual(values = col) +
  labs(y = 'Median Bray-Curtis dissimilarity', x = '', fill = '') +
  theme(legend.position = 'bottom') 
ggsave('out/metaphlan/mpa_boxplot_bray.png')  
