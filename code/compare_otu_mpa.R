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
library(microbiomics)

set.seed(96)
theme_set(theme_bw(base_size = 12) +
            theme(plot.title   = element_text(size = 12),
                  axis.title   = element_text(size = 12),
                  axis.text    = element_text(size = 12)))

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
sporulation_ability <- read.table('data/sporulation_ability/sporulation_ability2021.tsv', sep = '\t', header =T)
etoh_species <- read.table('data/longitudinal_shotgun/ethanol_resistant_species.tsv', sep = '\t', header =T)

bacteria <- read_metaphlan_table('~/projects/longitudinal_shotgun/data/metaphlan_abundance_table.txt', kingdom = "k__Bacteria", 
                                 lvl = 7, normalize = TRUE) %>% 
  rownames_to_column('name') %>% 
  pivot_longer(names_to = 'clade_name', values_to = 'value', cols = starts_with('k__')) %>% 
  mutate(name = str_remove_all(name, 'profiled_')) %>% 
  filter(name != 'MC013') %>% 
  left_join(metadata, by = join_by('name' == 'Group')) %>% 
  tidyr::separate_wider_delim(clade_name, delim = ".",
                              names = c('Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species')) %>%  
  mutate(Phylum = ifelse(Phylum == 'p__Firmicutes', 'p__Bacillota', Phylum), 
         Domain = str_remove_all(Domain, 'k__'), 
         Phylum = str_remove_all(Phylum, 'p__'), 
         Class = str_remove_all(Class, 'c__'), 
         Order = str_remove_all(Order, 'o__'), 
         Family = str_remove_all(Family, 'f__'), 
         Genus = str_remove_all(Genus, 'g__'), 
         Species = str_remove_all(Species, 's__'), 
         Species = str_replace_all(Species, '_', ' ')) %>% 
  left_join(select(sporulation_ability, PA, n_genes, sporulation_ability, Species), by = 'Species') %>% 
  left_join(select(etoh_species, Species, is_ethanol_resistant) %>%  
              mutate(Species = str_replace_all(Species, '_', ' ')), by = 'Species') %>% 
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
    TRUE ~ Phylum ))

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
col_phylum = c('#1F77B4', '#FF7F0E',  '#2CA02C',  '#D62728', '#9467BD', '#8C564B', '#f4d03f',  'pink', 'gray')


otu_long_all <-  readRDS('data/longitudinal_amplicons/otu_long_all.RDS') %>%
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
    TRUE ~ Phylum )) 

otu_rel <- otu_long_all %>% 
  group_by(Phylum, is_ethanol_resistant) %>%
  reframe(`16S amplicon data` = sum(rel_abund)/n_distinct(Group)*100) 

unique(otu_rel$Phylum)

mpa_rel <- bacteria %>%
  group_by(name, Phylum, is_ethanol_resistant) %>%
  reframe(total_abund = sum(value, na.rm = T)) %>% 
  group_by(Phylum, is_ethanol_resistant) %>% 
  reframe(`Metagenomic data` = sum(total_abund)/n_distinct(name)*100)
  
unique(mpa_rel$Phylum)

rel <- full_join(mpa_rel, otu_rel, by = c('Phylum', 'is_ethanol_resistant')) %>%
  pivot_longer(names_to = 'sample', values_to = 'rel_abund', cols = 3:4) %>% 
  mutate(phylum = ifelse(rel_abund < 0.1, '< 0.1 %', Phylum)) %>%  
  filter(!is.na(phylum))

rel %>%
  mutate(phylum = factor(phylum, levels = c('Bacillota', 'Bacteroidota', 'Actinomycetota', 'Pseudomonadota', 
                         'Verrucomicrobiota', 'Mycoplasmatota', 'Cyanobacteria', 'unclassified Bacteria', '< 0.1 %'))) %>% 
  ggplot(aes(y = is_ethanol_resistant, x = rel_abund, fill = phylum)) +
  geom_col() +
  scale_fill_manual(values = col_phylum) +
  facet_wrap(~sample, scales = 'free_x') +
  labs(x = 'Relative abundance [%]', y = '')
ggsave('out/compare_16s_meta/relative_abundance_col.png', dpi=400)

# Relative abudnance of these taxa in each sample
p_pyhlum <- otu_long_all %>%  
  filter(is_ethanol_resistant == c('Ethanol-resistant', 'Non ethanol-resistant')) %>% 
  filter(!Phylum %in% c('TM7', 'Deferribacteres', 'Lentisphaerota', 'Verrucomicrobiota')) %>% 
  group_by(Phylum, is_ethanol_resistant, person) %>% 
  reframe(rel_abund = mean(rel_abund, na.rm =T)) %>% 
  group_by(Phylum) %>%  
  reframe(test = list(wilcox.test(rel_abund[is_ethanol_resistant == "Ethanol-resistant"],
                                  rel_abund[is_ethanol_resistant == "Non ethanol-resistant"]))) %>% 
  mutate(p_value = sapply(test, function(x) x$p.value))
p_pyhlum

rel_otu_plot <- otu_long_all %>%  
  filter(is_ethanol_resistant == c('Ethanol-resistant', 'Non ethanol-resistant')) %>% 
  group_by(Phylum, is_ethanol_resistant, person) %>% 
  reframe(rel_abund = mean(rel_abund, na.rm =T)) %>% 
  ggplot(aes(y = Phylum, x = rel_abund, fill = is_ethanol_resistant)) +
  geom_boxplot()+
  geom_text(p_pyhlum, mapping = aes(y = Phylum, x = 1e-3, label = paste0('p =', signif(p_value, 2))), inherit.aes = F) + 
  scale_x_log10() +
  scale_fill_manual(values = c( '#f0a336', '#3CB371')) +
  labs(x = 'Relative abundance [%]', y = '', fill = '') +
  theme(legend.position = 'bottom', axis.text.y = element_text(face = 'italic')) 
rel_otu_plot
ggsave('out/compare_16s_meta/rel_otu_etoh_non.svg', dpi = 400)

# Metagenomic data 
p_phylum_mpa <- bacteria %>%  
  filter(is_ethanol_resistant == c('Ethanol-resistant', 'Non ethanol-resistant')) %>% 
  filter(!Phylum %in% c('Chloroflexota', 'Fusobacterium', 'Synergistota', 'unclassified Bacteria', 'Mycoplasmatota', 'Saccharibacteria')) %>% 
  group_by(Phylum, is_ethanol_resistant, person) %>% 
  reframe(value = mean(value, na.rm =T)) %>% 
  group_by(Phylum) %>%  
  reframe(test = list(wilcox.test(value[is_ethanol_resistant == "Ethanol-resistant"],
                                  value[is_ethanol_resistant == "Non ethanol-resistant"]))) %>% 
  mutate(p_value = sapply(test, function(x) x$p.value))
p_phylum_mpa

rel_mpa_plot <- bacteria %>% 
  filter(is_ethanol_resistant == c('Ethanol-resistant', 'Non ethanol-resistant')) %>% 
  group_by(Phylum, is_ethanol_resistant, person) %>% 
  reframe(value = mean(value, na.rm =T)) %>% 
  ggplot(aes(y = Phylum, x = value, fill = is_ethanol_resistant)) +
  geom_boxplot()+
  geom_text(p_phylum_mpa, mapping = aes(y = Phylum, x = 1e-6, label = paste0('p =', signif(p_value, 2))), inherit.aes = F) + 
  scale_x_log10() +
  scale_fill_manual(values = c( '#f0a336', '#3CB371')) +
  labs(x = 'Relative abundance [%]', y = '', fill = '') +
  theme(legend.position = 'bottom', axis.text.y = element_text(face = 'italic')) 
rel_mpa_plot
ggsave('out/compare_16s_meta/rel_species_etoh_non.svg', dpi = 400)

p1 <-ggarrange(rel_otu_plot + labs(title = '16S amplicon data'), 
          rel_mpa_plot + labs(title = 'Metagenomic data'), 
          common.legend = TRUE, legend = 'bottom')
p1
ggsave('out/compare_16s_meta/rel_both_etoh_non.svg', dpi = 400)

# For spore-forming and ethanol-resistant together! 
long_mpa <- readRDS('~/projects/longitudinal_amplicons/data/r_data/long_mpa.RDS')

rel_mpa_plot2 <- long_mpa %>%  
  filter(is_ethanol_resistant %in% c('Ethanol-resistant', 'Non ethanol-resistant'), 
         !is.na(sporulation_ability)) %>% 
  group_by(is_ethanol_resistant, sporulation_ability, person, Species) %>% 
  reframe(value = mean(value)) %>%  
  ggplot(aes(x = paste0(is_ethanol_resistant,' ' ,sporulation_ability), y = value)) +
  geom_boxplot() +
  scale_y_log10() +
  stat_compare_means(method = 'wilcox', p.adjust.method = "BH", comparisons = list(
    c("Ethanol-resistant Non-spore-former", "Non ethanol-resistant Non-spore-former"),
    c("Ethanol-resistant Non-spore-former", "Non ethanol-resistant Spore-former"),
    c("Ethanol-resistant Non-spore-former", "Ethanol-resistant Spore-former"),
    c("Ethanol-resistant Spore-former", "Non ethanol-resistant Spore-former"),
    c("Ethanol-resistant Spore-former", "Non ethanol-resistant Non-spore-former"),
    c("Non ethanol-resistant Spore-former", "Non ethanol-resistant Non-spore-former"))) +
  #scale_fill_manual(values = c( '#f0a336', '#3CB371', 'grey')) +
  #facet_wrap(~is_ethanol_resistant, scales = 'free') +
  labs(x = '', y = 'Relative abundance [%]') +
  theme_bw(base_size = 12) +
  theme(legend.position = 'bottom', 
        plot.title   = element_text(size = 12),
        axis.title   = element_text(size = 12),
        axis.text    = element_text(size = 11), 
        legend.text = element_text(size = 11)) 
rel_mpa_plot2

# OTUs simple
rel_otu_plot2 <- otu_long_all %>% 
  filter(is_ethanol_resistant %in% c('Ethanol-resistant', 'Non ethanol-resistant')) %>% 
  group_by(is_ethanol_resistant, person, name) %>% 
  reframe(rel_abund = mean(rel_abund)) %>%  
  ggplot(aes(x = is_ethanol_resistant, y = rel_abund)) +
  geom_boxplot() +
  scale_y_log10() +
  stat_compare_means(method = 'wilcox') +
  labs(x = '', y = 'Relative abundance [%]') +
  theme_bw(base_size = 12) +
  theme(legend.position = 'bottom', 
        plot.title   = element_text(size = 12),
        axis.title   = element_text(size = 12),
        axis.text    = element_text(size = 11), 
        legend.text = element_text(size = 11)) 


ggarrange(rel_otu_plot2 + labs(title = '16S amplicon data'), 
           rel_mpa_plot2 + labs(title = 'Metagenomic data'),
           align = 'v', widths = c(0.7, 1))
ggsave('out/compare_16s_meta/supplement_rel_both_simple.svg', dpi = 400)

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
alpha_long <- rbind(alpha_mpa %>%  mutate(sample = 'Metagenomic data'), alpha_otu %>% 
                      rename(name = Group, richness = S.obs) %>% 
                      select(name, richness, evenness, shannon) %>% 
                      mutate(sample = '16S amplicon data')) %>% 
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
  geom_point(size = 2) +
  facet_wrap(~person, scales = 'free') +
  labs(x='Day', y= 'Richness', color = 'Data', fill = 'Predefined\nevent', tag = 'A') +
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
  geom_point(size = 2) +
  facet_wrap(~person, scales = 'free') +
  labs(x='Day', y= 'Evenness', color = 'Data', fill = 'Predefined\nevent', tag = 'B') +
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

# Statistics if alpha diversity is in accrodance between sequencing methods
# Linear mixed-model - this takes into account the individuals (their points are not independant!)

alpha_lme <- alpha_otu %>% left_join(alpha_mpa, by = join_by('Group' == 'name')) 

library(lmerTest)
fit_rich <- lmer(S.obs ~ richness + day + (1 | person), data = alpha_lme) 
summary(fit_rich)

fit_even <- lmer(evenness.x ~ evenness.y + day + (1 | person), data = alpha_lme) 
summary(fit_even)


# Repeated measures correlation

library(rmcorr)
rmcorr(participant = as.factor(person), 
       measure1 = S.obs, 
       measure2 = richness, 
       data = alpha_lme)

rmcorr(participant = as.factor(person), 
       measure1 = evenness.x, 
       measure2 = evenness.y, 
       data = alpha_lme)




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
  geom_point(size = 3) +
  labs(color = 'Individual') +
  guides(color = guide_legend(nrow = 1)) +
  theme(legend.position = "bottom")
nmds_otu

# fit metadata 
meta1 <- metadata %>% 
  filter(Group %in% rownames(as.matrix(dist_otu))) %>% 
  column_to_rownames('Group') %>% 
  select(sample_type, age, height, weight, diet, food_supplements, 
         general_activity, digestion, household_members, living, environment, 
         diet14, antibiotics14, probiotics14, stress, bristol, 
         prebiotics14, medication14, moderate14, active14, dentist14, 
         vaccination14)

meta1[rownames(as.matrix(dist_otu)), ]
fit_otu <- envfit(dist_otu, meta1, permutations = 999, na.rm = TRUE)
fit_otu

# extract arrows
scores_env <- as.data.frame(scores(fit_otu, display = "vectors"))
scores_env$pval <- fit_otu$vectors$pvals

# plot
ord_otu_ar <- ggplot(beta_otu %>% filter(!is.na(person)), aes(x = NMDS1, y = NMDS2, color = person)) +
  geom_point(size = 3) +
  geom_segment(data = subset(scores_env, pval < 0.05),
               aes(x = 0, y = 0, xend = MA001, yend = MA002),
               arrow = arrow(length = unit(0.3, "cm")),
               color = "red") +
  geom_text(data = subset(scores_env, pval < 0.05),
            aes(x = MA001, y = MA002, label = rownames(scores_env %>% 
                                                         filter(pval < 0.05))),
            color = "black", vjust = -0.5) +
  labs(color = 'Individual') +
  theme(legend.position = "bottom") 
ord_otu_ar
ggsave('out/compare_16s_meta/OTUbeta_with_metadata.png')



# 
mpa_disp <- betadisper(dist_mpa, beta_mpa$person)
boxplot(mpa_disp)
ggsave('out/boxplot_betadisper_mpa.png')
anova(mpa_disp)
TukeyHSD(mpa_disp)

nmds_mpa <- beta_mpa %>% 
  filter(!is.na(person), biota == 'Bulk microbiota') %>% 
  ggplot(aes(x=NMDS1, y=NMDS2, color=person)) +
  geom_point(size = 3) +
  labs(color = 'Individual') +
  guides(color = guide_legend(nrow = 1)) +
  theme(legend.position = "bottom")
nmds_mpa

ggarrange(nmds_otu + labs(title ='16S amplicon data'), 
          nmds_mpa + labs(title ='Metagenomic data'), common.legend = T, legend = 'bottom') 
ggsave('out/compare_16s_meta/nmds_both.png', dpi = 600)

# fit metadata 
meta2 <- metadata %>% 
  filter(Group %in% rownames(as.matrix(dist_mpa))) %>% 
  column_to_rownames('Group') %>% 
  select(sample_type, age, height, weight, diet, food_supplements, 
         general_activity, digestion, household_members, living, environment, 
         diet14, antibiotics14, probiotics14, stress, bristol, 
         prebiotics14, medication14, moderate14, active14, dentist14, 
         vaccination14)

meta2[rownames(as.matrix(dist_mpa)), ]
fit_mpa <- envfit(dist_mpa, meta2, permutations = 999, na.rm = TRUE)
fit_mpa

# extract arrows
scores_env <- as.data.frame(scores(fit_mpa, display = "vectors"))
scores_env$pval <- fit_mpa$vectors$pvals

# plot
ord_mpa_ar <- ggplot(beta_mpa %>% filter(!is.na(person)), aes(x = NMDS1, y = NMDS2, color = person)) +
  geom_point(size = 3) +
  geom_segment(data = subset(scores_env, pval < 0.05),
               aes(x = 0, y = 0, xend = MA001, yend = MA002),
               arrow = arrow(length = unit(0.3, "cm")),
               color = "red") +
  geom_text(data = subset(scores_env, pval < 0.05),
            aes(x = MA001, y = MA002, label = rownames(scores_env %>% 
                                                         filter(pval < 0.05))),
            color = "black", vjust = -0.5) +
  labs(color = 'Individual') +
  theme(legend.position = "bottom") +
  guides(legend.guide = )
ord_mpa_ar
ggsave('out/compare_16s_meta/MPAbeta_with_metadata.png')


# What does time has to do with it? 
meta3 <- metadata %>% 
  filter(Group %in% rownames(as.matrix(dist_mpa))) %>% 
  column_to_rownames('Group') %>% 
  mutate(date = dmy(date), 
         season = ifelse(month(date) %in% c(12,1,2), "winter",
                         ifelse(month(date) %in% c(3,4,5), "spring",
                                ifelse(month(date) %in% c(6,7,8), "summer", "autumn")))) %>% 
  select(sample_type, season, age, height, weight, diet, food_supplements, 
         general_activity, digestion, household_members, living, environment, 
         diet14, antibiotics14, probiotics14, stress, bristol, 
         prebiotics14, medication14, moderate14, active14, dentist14, 
         vaccination14) 

fit_season <- envfit(dist_mpa, meta3[, c("season")], permutations = 999)
fit_season

# Both OTU and MPA
ggarrange(ord_otu_ar + labs(title ='16S amplicon data') + guides(color = guide_legend(nrow = 1)), 
          ord_mpa_ar + labs(title ='Metagenomic data') +  guides(color = guide_legend(nrow = 1)),
          common.legend = T, legend = 'bottom') 
ggsave('out/compare_16s_meta/beta_both_envfit.png', dpi = 400)


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

rel_both <- otu_long %>%  mutate(data = '16S amplicon data') %>% 
  rbind(mpa_long %>%  mutate(data = 'Metagenomic data')) %>% 
  mutate(biota = ifelse(biota == 'bulk microbiota', 'untreated sample', biota), 
         biota = factor(biota, levels = c('untreated sample', 'ethanol treated sample')), 
         Phylum = factor(Phylum, levels = c('Bacillota', 'Bacteroidota', 'Actinomycetota', 'Pseudomonadota',
                                            'Mycoplasmatota', 'Cyanobacteria', 'Verrucomicrobiota', 'unclassified Bacteria', '< 0.1%')))

ggplot(rel_both, aes(x = data, y = rel, fill = Phylum)) +
  geom_col() +
  facet_wrap(~biota, nrow = 1, scales = 'free_x') +
  scale_fill_manual(values = col_phylum) +
  labs(x = "", y = "Relative abundance [%]", fill = "Phylum") +
  theme(legend.text = element_text(face = 'italic'))
ggsave('out/compare_16s_meta/rel_abund.svg')
