library(readr)
library(tidyr)
library(dplyr)
library(tibble)
library(vegan)
library(lubridate)
library(ggpubr)
library(scales)

set.seed(96)
theme_set(theme_bw())

otutab <- readRDS('data/longitudinal_amplicons/otutab_ethanol_bulk.RDS')
taxtab <- readRDS('data/longitudinal_amplicons/taxtab.RDS')
metadata <- readRDS('data/longitudinal_amplicons/metadata.RDS')
norm_rel <- readRDS('data/longitudinal_amplicons/otutab_normrel.RDS') %>%
  left_join(select(metadata, Group, original_sample), by = 'Group')

long <- readRDS('data/longitudinal_amplicons/long_all.RDS')

# Colors to be used
col2 <- c('#3CB371', '#f0a336')

# Number of OTUs shared within person and between people
# OTU is present in a person only if it is present in more than 11 timepoints 
# Otherwise this analysis does not produce the same results, as too many OTUs are found by random 

n_otus <- long %>% 
  group_by(is_ethanol_resistant, person, Phylum, name) %>%
  reframe(present = sum(PA == 1)) %>%  
  filter(present > 11, Phylum %in% c('unclassified Bacteria', 'Pseudomonadota', 'Actinomycetota', 'Bacteroidota', 'Bacillota')) %>% 
  group_by(Phylum, is_ethanol_resistant) %>% 
  reframe(n = n_distinct(name))

n_otus_people <- long %>% mutate(time_point = as.integer(substr(Group, 3, 5))) %>%
  group_by(is_ethanol_resistant, person, Phylum, name) %>%
  reframe(present = sum(PA == 1)) %>%  
  filter(present >11, Phylum %in% c('unclassified Bacteria', 'Pseudomonadota', 'Actinomycetota', 'Bacteroidota', 'Bacillota')) %>% 
  mutate(present = ifelse(present > 0, 1, 0)) %>% 
  select(name, person, present, is_ethanol_resistant, Phylum) %>% 
  pivot_wider(names_from = 'person', values_from = 'present', values_fill = 0) %>%  
  mutate(n_people = A+B+C+D+E+F+G+H+I, 
         n_people = ifelse(n_people == 1, 'Present in a individual', 'Present in more than one individual')) %>% 
  group_by(is_ethanol_resistant, Phylum, n_people) %>% 
  reframe(n_otu = n_distinct(name)) %>% 
  left_join(n_otus, by = c('Phylum', 'is_ethanol_resistant')) %>% 
  mutate(per_otu = (n_otu/n)*100) 

n_otus_people$Phylum <- factor(n_otus_people$Phylum, levels = c('Bacillota', 'Bacteroidota',  
                                                          'Actinomycetota','Pseudomonadota', 
                                                          'unclassified Bacteria'))

n_otus_people %>% 
  ggplot(aes(x = is_ethanol_resistant, y = per_otu, fill = n_people))+
  geom_col() +
  scale_fill_manual(values = c('#F2933F', '#3F9EF2')) +
  facet_wrap(~Phylum, scales = 'free_y', nrow = 5) +
  labs(x = '', y = 'OTUs [%]', fill = '') +
  theme(legend.position = 'bottom') +
  guides(fill=guide_legend(nrow=2,byrow=TRUE))
ggsave('out/longitudinal_amplicons/present_one_more_people_phylum.png')

# not by Phylum 
n_otus <- long %>% 
  group_by(is_ethanol_resistant, person, name) %>%
  reframe(present = sum(PA == 1)) %>%  
  filter(present > 11) %>% 
  group_by(is_ethanol_resistant) %>% 
  reframe(n = n_distinct(name))

n_otus_people <- long %>% mutate(time_point = as.integer(substr(Group, 3, 5))) %>%
  group_by(is_ethanol_resistant, person, name) %>%
  reframe(present = sum(PA == 1)) %>%  
  filter(present > 11) %>% 
  mutate(present = ifelse(present > 0, 1, 0)) %>% 
  select(name, person, present, is_ethanol_resistant) %>% 
  pivot_wider(names_from = 'person', values_from = 'present', values_fill = 0) %>%  
  mutate(n_people = A+B+C+D+E+F+G+H+I, 
         n_people = ifelse(n_people == 1, 'Present in a individual', 'Present in more than one individual')) %>% 
  group_by(is_ethanol_resistant, n_people) %>% 
  reframe(n_otu = n_distinct(name)) %>% 
  left_join(n_otus, by = 'is_ethanol_resistant') %>% 
  mutate(per_otu = (n_otu/n)*100) 

n_otus_people %>% 
  ggplot(aes(x = is_ethanol_resistant, y = per_otu, fill = n_people))+
  geom_col() +
  scale_fill_manual(values = c('#F2933F', '#3F9EF2')) +
  labs(x = '', y = 'OTUs [%]', fill = '') +
  theme(legend.position = 'bottom') +
  guides(fill=guide_legend(nrow=2,byrow=TRUE))
ggsave('out/longitudinal_amplicons/present_one_more_people.png')

  


