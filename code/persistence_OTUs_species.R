# Persistence of OTUs and species within individual 

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

long_otu <- readRDS('data/longitudinal_amplicons/long_all.RDS')

persistence_otu <- long_otu %>%
  mutate(time_point = as.integer(substr(Group, 3, 5))) %>%
  group_by(is_ethanol_resistant, person, name) %>%
  reframe(all_timepoints = n(), 
          timepoints_present = sum(PA == 1)) %>%
  # Remove 'singletons'
  filter(timepoints_present > 1) %>%
  mutate(prevalence = (timepoints_present/all_timepoints)*100) %>%
  group_by(person, is_ethanol_resistant) %>%
  mutate(no_otus = n_distinct(name)) %>%
  ungroup() %>%
  group_by(is_ethanol_resistant, person, prevalence) %>%
  reframe(no_otus2 = n_distinct(name), 
          per_otus = (no_otus2/no_otus) *100) %>%
  mutate(sporulation_ability = 'Non-spore-former') %>% 
  unique()
persistence_otu

# Shotgun 
long_mpa <- readRDS('data/longitudinal_shotgun/long_mpa.RDS') %>% 
  filter(is_ethanol_resistant %in% c('Ethanol-resistant', 'Non ethanol-resistant'), 
         biota == 'untreated sample') 

persistence_mpa <- long_mpa %>%
  mutate(time_point = as.integer(substr(name, 3, 5)), 
         pa = ifelse(value > 0, 1, 0)) %>%
  group_by(is_ethanol_resistant, sporulation_ability, person, Species) %>%
  reframe(all_timepoints = n(), 
          timepoints_present = sum(pa)) %>%
  # Remove 'singletons'
  filter(timepoints_present > 1) %>%
  mutate(prevalence = (timepoints_present/all_timepoints)*100) %>%
  group_by(person, is_ethanol_resistant) %>%
  mutate(no_otus = n_distinct(Species)) %>%
  ungroup() %>%
  group_by(is_ethanol_resistant, sporulation_ability, person, prevalence) %>%
  reframe(no_otus2 = n_distinct(Species), 
          per_otus = (no_otus2/no_otus) *100) %>%
  unique()

persistence <- rbind(persistence_mpa %>%  mutate(method = 'metagenomic data'), 
                     persistence_otu %>%  mutate(method = '16S amplicon data')) %>% 
  left_join(select(long_mpa, is_ethanol_resistant, sporulation_ability) %>% unique(), by = c('is_ethanol_resistant', 'sporulation_ability'), 
            relationship = "many-to-many") 

# Original plot 
plot_persist_mpa <- ggplot(persistence_mpa, aes(x = prevalence, y = per_otus, color = is_ethanol_resistant, linetype = sporulation_ability)) +
  #geom_point(size = 3, alpha = 0.3) +
  geom_smooth(method = "loess", formula = y ~ x, se = F, linewidth = 1.5) +
  scale_color_manual(values = c('#f0a336', '#3CB371'))+
  labs(x='Within-individual persistence\n [% of time points present]', y= 'Taxa per individual\n [% of taxa]', color = '', linetype = '',  title = 'metagenomic data') +
  guides(color = guide_legend(nrow = 2), 
         linetype = guide_legend(nrow = 2)) +
  theme(legend.position = 'bottom') 
plot_persist_mpa

plot_persist_otu <- ggplot(persistence_otu, aes(x = prevalence, y = per_otus, color = is_ethanol_resistant, linetype = sporulation_ability)) +
  #geom_point(size = 3, alpha = 0.3) +
  geom_smooth(method = "loess", formula = y ~ x, se = F, linewidth = 1.5) +
  scale_color_manual(values = c('#f0a336', '#3CB371'))+
  labs(x='Within-individual persistence\n [% of time points present]', y= 'Taxa per individual\n [% of taxa]', color = '', linetype = '', title = '16S amplicon data') +
  guides(color = guide_legend(nrow = 2)) +
  theme(legend.position = 'bottom') 
plot_persist_otu

# Alternative plot (not line)
plot_persist_mpa2 <- persistence_mpa %>%  
  filter(!is.na(sporulation_ability)) %>% 
  ggplot(aes(x = prevalence, y = per_otus, color = is_ethanol_resistant)) +
  geom_point(size = 4, alpha = .8) +
  geom_smooth(method = "loess", formula = y ~ x, se = F, linewidth = 1, alpha = .4) +
  scale_color_manual(values = c('#f0a336', '#3CB371'))+
  facet_wrap(~sporulation_ability) +
  labs(x='Within-individual persistence\n [% of time points present]', y= 'Species per individual\n [% of species]', 
       color = '', title = 'Metagenomic data') +
  guides(color = guide_legend(nrow = 1)) +
  theme(legend.position = 'bottom') 
plot_persist_mpa2

plot_persist_otu2 <- ggplot(persistence_otu, aes(x = prevalence, y = per_otus, color = is_ethanol_resistant)) +
  geom_point(size = 4, alpha = .8) +
  geom_smooth(method = "loess", formula = y ~ x, se = F, linewidth = 1, alpha = .4) +
  scale_color_manual(values = c('#f0a336', '#3CB371'))+
  labs(x='Within-individual persistence\n [% of time points present]', y= 'OTUs per individual\n [% of OTUs]', color = '', title = '16S amplicon data') +
  guides(color = guide_legend(nrow = 1)) +
  theme(legend.position = 'bottom') 
plot_persist_otu2

ggarrange(plot_persist_otu2, plot_persist_mpa2,
          common.legend = TRUE, legend = 'bottom',
          widths = c(.8, 1))
ggsave('out/compare_16s_meta/persistence_alternative.svg', dpi=400)

# Statisitcs, is there really more ethanol resistant OTUs present at more time-points than ethanol non-resistant 

# Linear mixed-effects model does ethanol-resistency and sporulation influence persistence within a person? 
# This is testing if mean is different.. 

# For OTUs 
persistence_otu_stat <- long_otu %>% 
  mutate(time_point = as.integer(substr(Group, 3, 5))) %>%
  group_by(is_ethanol_resistant, person, name) %>%
  reframe(all_timepoints = n(), 
          timepoints_present = sum(PA)) %>%
  # Remove 'singletons'
  filter(timepoints_present > 1) %>%
  mutate(prevalence = (timepoints_present/all_timepoints)*100) 

wilcox.test(prevalence ~ is_ethanol_resistant, data = persistence_otu_stat,
            exact = FALSE, conf.int = TRUE) 
# W = 12883706, p-value < 2.2e-16

# Metagenomic data 
persistence_mpa_stat <- long_mpa %>%
  filter(biota == 'untreated sample', !is.na(sporulation_ability)) %>% 
  mutate(time_point = as.integer(substr(name, 3, 5)), 
         pa = ifelse(value > 0, 1, 0)) %>%
  group_by(is_ethanol_resistant, sporulation_ability, person, Species) %>%
  reframe(all_timepoints = n(), 
          timepoints_present = sum(pa)) %>%
  # Remove 'singletons'
  filter(timepoints_present > 1) %>%
  mutate(prevalence = (timepoints_present/all_timepoints)*100) %>% 
  mutate(group = interaction(is_ethanol_resistant, sporulation_ability, sep = " "))

pairwise.wilcox.test(x = persistence_mpa_stat$prevalence,
                     g = persistence_mpa_stat$group,
                     p.adjust.method = "BH" )

# Pairwise comparisons using Wilcoxon rank sum test with continuity correction 
# 
# data:  persistence_mpa_stat$frac_high_persist and persistence_mpa_stat$group 
# 
#                                        Ethanol-resistant Non-spore-former
# Non ethanol-resistant Non-spore-former 0.0408                            
# Ethanol-resistant Spore-former         0.0032                            
# Non ethanol-resistant Spore-former     0.0622                            
#                                         Non ethanol-resistant Non-spore-former
# Non ethanol-resistant Non-spore-former  -                                     
# Ethanol-resistant Spore-former          0.0040                                
# Non ethanol-resistant Spore-former      0.5074                                
#                                        Ethanol-resistant Spore-former
# Non ethanol-resistant Non-spore-former  -                             
# Ethanol-resistant Spore-former          -                             
# Non ethanol-resistant Spore-former      0.0024                 

# 
# P value adjustment method: BH 