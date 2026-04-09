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
metadata <- read.table('data/metadata.csv', sep = ';', header = T) %>% 
  mutate(biota = ifelse(biota == 'bulk microbiota', 'untreated sample', biota))

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
  filter(!is.na(biota)) %>% 
  mutate(biota = ifelse(biota == 'bulk microbiota', 'untreated sample', biota), 
         biota = factor(biota, levels =c('untreated sample', 'ethanol treated sample'))) 


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
  geom_point(size = 3) +
  scale_color_manual(values = c('#3CB371', '#f0a336')) +
  facet_wrap(~biota, scales = 'free') +
  geom_abline()+  
  #scale_x_continuous(limits = c(50, 200)) +
  #scale_y_continuous(limits = c(50, 200)) +
  labs(x = '# genera / sample\n [amplicon data]', y = '# genera / sample\n [shotgun metagenomic data]',  color = '') +
  theme(legend.position = 'none')
ggsave('out/compare_genera_original.png', dpi = 600)

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

alpha_mpa <- readRDS('data/longitudinal_shotgun/alpha_diveristy.RDS') %>% select(name, richness, shannon, evenness)
alpha_otu <- readRDS('data/longitudinal_amplicons/alpha.RDS')

alpha <- full_join(alpha_otu, alpha_mpa, by = join_by('Group' == 'name'))

events <- read.table('data/extreme_event_data.csv', sep = ',', header = T)

alpha %>% filter(biota == 'Bulk microbiota') %>% 
  pivot_longer(names_to = 'sample', values_to = 'shannon', cols = starts_with('shannon')) %>% 
  mutate(sample = ifelse(sample == 'shannon.x', 'Amplicon data', 'Shotgun metagenomic data')) %>% 
  ggplot(aes(x=day, y=shannon, color = sample)) +
  geom_rect(data = events, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = extremevent_type), inherit.aes = FALSE,
            alpha = 0.6) +
  scale_fill_manual(values = c('white','#d94343', '#d98e43', '#f1f011', '#0c9910', '#3472b7', '#7934b7', '#b73485', '#0f5618')) +
  scale_color_manual(values = c('skyblue', 'violet')) +
  geom_line(linewidth=1.5) +
  # geom_line(data=alpha_meta %>% dplyr::select(-person) %>% filter(biota == 'Ethanol treated sample'), 
  #           aes(group=person2), color= cole, linewidth=0.5, alpha=0.5)+
  # geom_line(data=alpha_meta %>% filter(biota == 'Bulk microbiota'),
  #           aes(color=person), color= colm, linewidth=1.2) +
  # geom_line(data=alpha_meta %>% filter(biota == 'Ethanol treated sample'), 
  #           color=cole, linewidth=1.2) +
  facet_wrap(~person, scales = 'free_x') +
  labs(x='Day', y= 'Shannon', color = 'Data type', fill = 'Event') +
  theme(legend.position = 'bottom')
ggsave('out/alpha_otu_mpa.svg', dpi = 600)


# Richness and evenness
alpha_long <- rbind(alpha_mpa %>%  mutate(sample = 'Metagenomics'), alpha_otu %>% 
                      rename(name = Group, richness = S.obs) %>% 
                      select(name, richness, evenness, shannon) %>% 
                      mutate(sample = '16 AS')) %>% 
  left_join(metadata, by = join_by('name' == 'Group'))

# Correlation by person? 
richness_plot <- alpha_long %>%  
  filter(biota == 'untreated sample') %>% 
  ggplot(aes(x=day, y=richness, color = sample)) +
  geom_rect(data = events, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = extremevent_type), inherit.aes = FALSE,
            alpha = 0.6) +
  scale_fill_manual(values = c('white','#d94343', '#d98e43', '#f1f011', '#0c9910', '#3472b7', '#7934b7', '#b73485', '#0f5618')) +
  scale_color_manual(values = c('skyblue', 'violet')) +
  geom_line(linewidth=1.5) +
  facet_wrap(~person, scales = 'free') +
  labs(x='Day', y= 'Richness', color = 'Data', fill = 'Event', title = 'Richness') +
  theme(legend.position = 'bottom')
richness_plot

evenness_plot <- alpha_long %>%  
  filter(biota == 'untreated sample') %>% 
  ggplot(aes(x=day, y=evenness, color = sample)) +
  geom_rect(data = events, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = extremevent_type), inherit.aes = FALSE,
            alpha = 0.6) +
  scale_fill_manual(values = c('white','#d94343', '#d98e43', '#f1f011', '#0c9910', '#3472b7', '#7934b7', '#b73485', '#0f5618')) +
  scale_color_manual(values = c('skyblue', 'violet')) +
  geom_line(linewidth=1.5) +
  facet_wrap(~person, scales = 'free') +
  labs(x='Day', y= 'Evenness', color = 'Data', fill = 'Event', title = 'Evenness') +
  theme(legend.position = 'bottom')
evenness_plot

ggarrange(richness_plot, evenness_plot, 
          common.legend = TRUE, legend = 'bottom', 
          ncol = 1)
ggsave('out/compare_16s_meta/richness_evenness.svg')

# 
ggarrange(alpha %>% filter(!is.na(biota)) %>% 
            ggplot(aes(y = S.obs, x = richness,  color = biota)) +
            geom_point(size = 3) +
            geom_abline() +
            scale_color_manual(values = c('#3CB371', '#f0a336')) +
            facet_wrap(~biota, scales = 'free') +
            labs(y = 'Richness [amplicon data]', x = 'Richness [shotgun data]', tag = 'A') +
            theme(legend.position = 'none'), 
          alpha %>% filter(!is.na(biota)) %>% 
            ggplot(aes(y = evenness.x, x = evenness.y, color = biota)) +
            geom_point(size = 3) +
            geom_abline() +
            scale_color_manual(values = c('#3CB371', '#f0a336')) +
            facet_wrap(~biota, scales = 'free') +
            labs(y = 'Evenness [amplicon data]', x = 'Evenness [shotgun data]', tag = 'B') +
            theme(legend.position = 'none'), nrow = 2)
ggsave('out/richness_evenntess_compare.png') 


# Plot for presentation 
ggarrange(alpha %>% filter(biota == 'Bulk microbiota') %>% 
            ggplot(aes(y = shannon.x, x = shannon.y)) +
            geom_point(size = 3, aes(color = person)) +
            geom_smooth(method = 'lm') +
            stat_cor() +
            labs(y = 'Shannon [16S amplicon data]', x = 'Shannon [metagenomic data]', title = 'Shannon') +
            theme(legend.position = 'none'), 
          alpha %>% filter(biota == 'Bulk microbiota') %>% 
            ggplot(aes(y = evenness.x, x = evenness.y)) +
            geom_point(size = 3, aes(color = person)) +
            geom_smooth(method = 'lm') +
            stat_cor() +
            labs(y = 'Evenness [16S amplicon data]', x = 'Evenness [metagenomic data]', title = 'Evenness') +
            theme(legend.position = 'none'), nrow = 2)

ggsave('out/compare_16s_meta/richness_evenntess_compare.png') 

alpha %>% filter(biota == 'Bulk microbiota') %>% 
  ggplot(aes(y = shannon.x, x = shannon.y)) +
  geom_point(size = 3, aes(color = person)) +
  geom_smooth(method = 'lm') +
  stat_cor() +
  labs(y = 'Shannon [16S amplicon data]', x = 'Shannon [metagenomic data]') +
  theme(legend.position = 'none')
ggsave('out/compare_16s_meta/shannon_corr.png')

# Beta diveristy 
beta_mpa <- readRDS('data/longitudinal_shotgun/nmds_mpa_positions.rds')
beta_otu <- readRDS('data/longitudinal_amplicons/nmds_otu_positions.RDS')

dist_otu <- readRDS('data/longitudinal_amplicons/dist_otu.RDS')
dist_mpa <- readRDS('data/longitudinal_shotgun/dist_mpa.RDS')

# Betadisper 
otu_disp <- betadisper(dist_otu, beta_otu$person)
boxplot(otu_disp)
ggsave('out/boxplot_betadisper_otu.png')
anova(otu_disp)
TukeyHSD(otu_disp)

nmds_otu <- beta_otu %>%
  ggplot(aes(x=NMDS1, y=NMDS2, color=person)) +
  geom_point(size = 5) +
  labs(color = 'Individual') +
  guides(color = guide_legend(nrow = 1)) +
  theme(legend.position = "bottom")
nmds_otu

# 
mpa_disp <- betadisper(dist_mpa, beta_mpa$person)
boxplot(mpa_disp)
ggsave('out/boxplot_betadisper_mpa.png')
anova(mpa_disp)
TukeyHSD(mpa_disp)


nmds_mpa <- beta_mpa %>% 
  filter(!is.na(person), biota == 'Bulk microbiota') %>% 
  ggplot(aes(x=NMDS1, y=NMDS2, color=person)) +
  geom_point(size = 5) +
  labs(color = 'Individual') +
  guides(color = guide_legend(nrow = 1)) +
  theme(legend.position = "bottom")
nmds_mpa

ggarrange(nmds_otu + labs(tag ='A'), 
          nmds_mpa + labs(tag ='B'), common.legend = T, legend = 'bottom') 
ggsave('out/nmds_both.png')

ggarrange(nmds_otu + labs(title ='16S amplicon sequencing'), 
          nmds_mpa + labs(title ='Metagenomic sequencing'), common.legend = T, legend = 'bottom') 
ggsave('out/compare_16s_meta/nmds_both_presentation.svg')

# Relative abundance of bacteria for OTUs and MPA 
otutab <- readRDS('data/longitudinal_amplicons/otutab_ethanol_bulk.RDS')
taxtab <- readRDS('data/longitudinal_amplicons/taxtab.RDS')
col_phylum = c('#1F77B4', '#FF7F0E',  '#2CA02C',  '#D62728', '#9467BD', '#8C564B', '#f4d03f', '#0f5618', '#3127F5')


otu_long <- as.data.frame(otutab) %>% 
  rownames_to_column('Group') %>%  
  pivot_longer(names_to = 'name', values_to = 'value', cols = starts_with('Otu')) %>% 
  left_join(metadata, by = 'Group') %>% 
  left_join(taxtab, by = 'name') %>% 
  group_by(biota) %>% 
  mutate(rel_abund = value/sum(value) *100) %>% 
  group_by(biota, Phylum) %>% 
  reframe(rel = sum(rel_abund)) %>% 
  mutate(Phylum = ifelse(rel < 0.1, '< 0.1%', Phylum), 
         Phylum = ifelse(Phylum == 'Tenericutes', 'Mycoplasmatota', Phylum)) %>%  
  group_by(biota, Phylum) %>% 
  reframe(rel = sum(rel))

otu_long$Phylum <- factor(otu_long$Phylum, levels = c('Bacillota', 'Bacteroidota', 'Actinomycetota', 'Pseudomonadota',
                                                    'Verrucomicrobiota', 'Mycoplasmatota' ,'unclassified Bacteria', '< 0.1%'))
# Relative abundance plots for comparison
rel_otu <- otu_long %>%
  ggplot(aes(x = biota, y = rel, fill = Phylum)) +
  geom_col() +
  scale_fill_manual(values = col_phylum) +
  labs(x = '', y= 'Relative abundance [%]', fill = 'Phylum')
rel_otu 

#mpa 
mpa_long <- bacteria %>% 
  group_by(name, biota, Phylum) %>%  
  reframe(rel = sum(value)) %>% 
  group_by(biota, Phylum) %>%  
  reframe(rel = mean(rel)) %>% 
  mutate(Phylum = ifelse(rel < 0.1, '< 0.1%', Phylum)) %>%  
  mutate(Phylum = case_when(
        Phylum == 'Actinobacteria' ~ 'Actinomycetota',
        Phylum == 'Candidatus_Saccharibacteria' ~ 'Saccharibacteria',
        Phylum == 'Chloroflexi' ~ 'Chloroflexota',
        Phylum == 'Proteobacteria' ~ 'Pseudomonadota',
        Phylum == 'Lentisphaerae' ~ 'Lentisphaerota',
        Phylum == 'Fusobacteria' ~ 'Fusobacterium',
        Phylum == 'Verrucomicrobia' ~ 'Verrucomicrobiota',
        Phylum == 'Bacteria_unclassified' ~ 'unclassified Bacteria',
        Phylum == 'Candidatus_Melainabacteria' ~ 'Cyanobacteria',
        Phylum == 'Synergistetes' ~ 'Synergistota',
        Phylum == 'Tenericutes' ~ 'Mycoplasmatota',
        TRUE ~ Phylum )) %>%  
  group_by(biota, Phylum) %>% 
  reframe(rel = sum(rel))
  

mpa_long$Phylum <- factor(mpa_long$Phylum, levels = c('Bacillota', 'Bacteroidota', 'Actinomycetota', 'Pseudomonadota',
                                                      'Mycoplasmatota', 'Cyanobacteria', 'Verrucomicrobiota','< 0.1%'))

rel_both <- otu_long %>%  mutate(data = 'Amplicon data') %>% 
  rbind(mpa_long %>%  mutate(data = 'Shotgun data')) %>% 
  mutate(biota = ifelse(biota == 'bulk microbiota', 'untreated sample', biota), 
         biota = factor(biota, levels = c('untreated sample', 'ethanol treated sample')))

ggplot(rel_both, aes(x = data, y = rel, fill = Phylum)) +
  geom_col() +
  facet_wrap(~biota, nrow = 1, scales = 'free_x') +
  scale_fill_manual(values = col_phylum) +
  labs(x = "", y = "Relative abundance [%]", fill = "Phylum") 
ggsave('out/rel_abund.png')
