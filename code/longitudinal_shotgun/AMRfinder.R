library(dplyr)
library(readr)
library(tibble)
library(tidyr)
library(ggplot2)
library(vegan)
library(stringr)
library(forcats)
library(ggpubr)
library(lubridate)

set.seed(96)
theme_set(theme_bw())

metadata <- read_csv2('~/projects/longitudinal_shotgun/data/metadata.csv') 

# Read in ARMfinder data 
amr_pre <- read_tsv('~/projects/longitudinal_shotgun/data/AMR/AMRfinder/amrfinder_results.tsv') %>%
  rename('contig_id' = `Contig id`, 'ARG' = `Element symbol`) %>%
  filter(`% Coverage of reference` > 80 & `% Identity to reference` > 50) %>%  
  select(contig_id, ARG, Class)

unique(amr_pre$ARG)

mechanisms <- read_csv2('~/projects/longitudinal_shotgun/data/AMR/AMRfinder/mechanism_resistence.csv')

# Run this lines only the first time
# Extract names of contigs for which ARGs were found
# orf_ids = amr_pre$contig_id
# writeLines(orf_ids, 'data/intermediate/AMRfinder_ORF.list')

# Before you run this script forward, make sure to downsize the 10.<>.mapcount file with Python script code/filter_AMRfinder_ORFs.sh
# bash filter_AMRfinder_ORFs.sh
contig_amr_coverage <- read.table('~/projects/longitudinal_shotgun/data/sqm_tables/19.filtered_AMRfinder.contigtable', header = TRUE, sep = '\t') %>%
  rename('contig_id' = 'Contig.ID') %>%
  select(-Bin.ID, -starts_with('Raw.'), -starts_with('TPM')) %>%
  pivot_longer(names_to = 'Group', values_to = 'coverage', cols = starts_with('Cov')) %>%
  #pivot_longer(names_to = 'Group2', values_to  = 'TPM', cols = starts_with('TPM')) %>%
  mutate(Group = str_sub(Group, 10, 14))

contig_amr_TPM <- read.table('~/projects/longitudinal_shotgun/data/sqm_tables/19.filtered_AMRfinder.contigtable', header = TRUE, sep = '\t') %>%
  rename('contig_id' = 'Contig.ID') %>%
  select(-Bin.ID, -starts_with('Raw.'), -starts_with('Cov')) %>%
  pivot_longer(names_to = 'Group', values_to  = 'TPM', cols = starts_with('TPM')) %>%
  mutate(Group = str_sub(Group, 5, 9)) 


amr <- inner_join(contig_amr_coverage, contig_amr_TPM, by = c('contig_id', 'Tax', 'Group', 'Disparity', 'GC.perc', 'Length', 'Num.genes')) %>% 
  full_join(amr_pre, by = 'contig_id', relationship = 'many-to-many') %>% 
  left_join(metadata, by ='Group', relationship = "many-to-many") %>%
  left_join(mechanisms, by = 'ARG') %>% 
  filter(TPM > 0 | biota == 'bulk microbiota') %>% 
  filter(!is.na(Class))
   

amr %>% summarise(n_distinct(ARG))

amr %>% summarise(n_distinct(Class))

amr_simplyfied <- amr %>%
  group_by(Class) %>%
  mutate(sum_unique = n_distinct(ARG)) %>% 
  ungroup() %>% 
  filter(sum_unique > 5) %>% 
  select(-sum_unique)

amr_simplyfied$mechanism_resistence <- factor(amr_simplyfied$mechanism_resistence, 
                                              levels = c('antibiotic efflux', 
                                                          'antibiotic inactivation', 
                                                          'antibiotic target protection', 
                                                          'antibiotic target alteration', 
                                                          'antibiotic target replacement'))
amr_simplyfied %>% summarise(n_distinct(ARG))

amr_simplyfied %>% summarise(n_distinct(Class))

# unique classes 
unique <- amr %>%
  group_by(Class) %>%
  summarise(sum_unique = n_distinct(ARG)) %>%
  ggplot(aes(x = sum_unique, y = reorder(Class, sum_unique))) +
  geom_col() +
  labs(x = '# unique genes', y = 'Class of ARG')
unique
ggsave('out/ARGs/AMRf_class_unique.png')

# tpm of each unique class 
tpm <- amr %>%
  group_by(Class) %>%
  reframe(sum_TPM = sum(TPM),.groups = 'drop') %>%
  ggplot(aes(x = sum_TPM, y = reorder(Class, sum_TPM))) +
  geom_col() +
  scale_x_log10() +
  labs(y = 'Class of ARG', x = 'log(TPM)') 
  
  #caption = 'TPM = a feature (be it a transcript, a gene or a functional category) 
  #     the number of times that we would find that feature when randomly sampling 1 million features, 
  #     given the abundances of the different features in our sample')
tpm
ggsave('out/ARGs/AMRf_TPM_unique_class_ARGs_log.png')

ggarrange(unique, tpm, common.legend = TRUE)
ggsave('out/ARGs/AMRf_TPM_class_ARGs.png')

# mechanism of resistence 
amr %>%
  group_by(mechanism_resistence) %>%
  reframe(sum_TPM = sum(TPM),.groups = 'drop') %>%
  ggplot(aes(x = sum_TPM, y = reorder(mechanism_resistence, sum_TPM))) +
  geom_col() +
  labs(y = 'Mechanism of resistence', x = 'TPM') 
ggsave('out/ARGs/AMRf_mechanisms_tpm.png')

amr %>%
  group_by(mechanism_resistence) %>%
  reframe(unique = n_distinct(ARG),.groups = 'drop') %>%
  ggplot(aes(x = unique, y = reorder(mechanism_resistence, unique))) +
  geom_col() +
  scale_x_log10() +
  labs(y = 'Mechanism of resistence', x = '# unique ARGs') 
ggsave('out/ARGs/AMRf_mechanisms_unique.png')

# Number of unique mechanisms colored by class
amr_tpm <- amr %>%
  group_by(mechanism_resistence, Class) %>%
  reframe(sum_TPM = round(sum(TPM))) %>% 
  ggplot(aes(x = sum_TPM, y = reorder(Class, sum_TPM), fill = mechanism_resistence)) +
  geom_col() +
  geom_label(aes(x = sum_TPM, label = sum_TPM), size = 3) +
  labs(y = '', x = 'TPM', fill = 'Mechanism of resistence') +
  theme(legend.position = 'bottom') +
  guides(fill = guide_legend(nrow = 1))
amr_tpm

amr_unique <- amr %>%
  group_by(mechanism_resistence, Class) %>%
  reframe(unique = n_distinct(ARG)) %>%
  group_by(Class) %>% 
  mutate(sum_unique = sum(unique)) %>% 
  ggplot(aes(x = unique, y = reorder(Class, sum_unique), fill = mechanism_resistence)) +
  geom_col() +
  geom_label(aes(x = unique, label = unique), size = 3) +
  labs(y = '', x = '# unique ARGs', fill = 'Mechanism of resistence') +
  theme(legend.position = 'bottom') +
  guides(fill = guide_legend(nrow = 1))
amr_unique

ggarrange(amr_unique + labs(tag = 'A'), 
          amr_tpm + labs(tag = 'B'), 
          common.legend = T, legend = 'bottom', 
          nrow = 1)

ggsave('out/ARGs/unique_tpm_args_mechanism.png', dpi = 600)
library(svglite)
ggsave('out/ARGs/unique_tpm_args_mechanism.svg', dpi = 600)


# TPM of different ARG classes through time 
# Taking into account events 
event_data <- read.table('data/extreme_event_data.csv', sep = ',', header = TRUE) 

amr_simplyfied %>%
  group_by(person, day, Class) %>%
  reframe(sum_TPM = sum(TPM), extremevent_type) %>%
  ggplot(aes(x = day, y = sum_TPM)) +
  geom_rect(data = event_data, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = extremevent_type), inherit.aes = FALSE,
            alpha = 0.6) +
  scale_fill_manual(values = c('white','#d94343', '#d98e43', '#f1f011', '#0c9910', '#3472b7', '#7934b7', '#b73485', '#0f5618')) +
  geom_point(size = 2, aes(color  = Class)) +
  geom_line(linewidth = 1, aes(color  = Class)) +
  facet_wrap(~person, scales = 'free_y') +
  labs(y = 'TPM [log10]', x = 'Day', color = 'Class of ARG', fill = 'Event type') 
ggsave('out/ARGs/AMRf_extreme_event_TPM.png')

# Number of different ARGs per individual through time 
amr_simplyfied %>% filter(TPM > 0) %>%
  group_by(person, day, Class) %>%
  summarise(sum_unique = n_distinct(ARG), .groups = 'drop') %>%
  filter(sum_unique > 0) %>%
  ggplot(aes(x = day, y = sum_unique)) +
  geom_rect(data = event_data, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = extremevent_type), inherit.aes = FALSE,
            alpha = 0.6) +
  scale_fill_manual(values = c('white','#d94343', '#d98e43', '#f1f011', '#0c9910', '#3472b7', '#7934b7', '#b73485', '#0f5618')) +
  geom_point(size = 2, aes(color = Class)) +
  geom_line(linewidth=1, aes(color = Class)) +
  facet_wrap(~person, scales = 'free_y') +
  labs(y = '# unique genes', x = 'Day', color = 'Class of ARG', fill = 'Event type') 
ggsave('out/ARGs/AMRf_unique_class_ARG.png')

# What is the taxonomy of different ARG classes 
amr_tax <- amr_simplyfied %>%
  mutate(Phylum = str_extract(Tax, "p_[^;]+"),
         Phylum = substr(Phylum, 3, 25), 
         Kingdom = str_extract(Tax, 'k_[^;]+'), 
         Class2 = str_extract(Tax, 'c_[^;]+'))

# How much reads mapped from each phyla? 
amr_tax %>%
  group_by(Class, Phylum) %>% 
  reframe(TPM = sum(TPM)) %>% 
  ggplot(aes (x = TPM, y = reorder(Class, TPM), fill = Phylum)) +
  geom_col() +
  #geom_label(aes(label = TPM)) +
  scale_fill_manual(values = c('#d94343', '#0c9910','#3472b7', '#b73485', '#f1f011', 'lightgrey' )) +
  theme(legend.position = 'bottom') +
  labs(y = '')
ggsave('out/ARGs/AMRf_taxonomy_class_TPM_simple.png')

# Which phylum has the highest number of unique ARGs? 
amr_tax$Phylum <- factor(amr_tax$Phylum, levels = c('NA', 'Pseudomonadota', 'Actinomycetota', 'Bacteroidota', 'Bacillota'))

unique_tax <- amr_tax %>%
  group_by(Class, Phylum) %>%  
  reframe(unique = n_distinct(ARG)) %>%  
  ggplot(aes (x = unique, y = Phylum, fill = Class)) +
  geom_col() +
  #scale_fill_manual(values = c('#d94343', '#0c9910','#3472b7', '#b73485', '#f1f011', 'lightgrey' )) +
  theme(legend.position = 'bottom') +
  labs(y = '', x = '# unique ARGs')
ggsave('out/ARGs/AMRf_taxonomy_class_unique_simple.png')

# Which phylum has the highest abundance of ARGs? 
tpm_tax <- amr_tax %>%
  group_by(Class, Phylum) %>%  
  reframe(TPM = sum(TPM)) %>%  
  ggplot(aes (x = TPM, y = Phylum, fill = Class)) +
  geom_col() +
  #scale_fill_manual(values = c('#d94343', '#0c9910','#3472b7', '#b73485', '#f1f011', 'lightgrey' )) +
  theme(legend.position = 'bottom', axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  labs(y = '')
ggsave('out/ARGs/AMRf_taxonomy_class_TPM_simple_unique.png')

ggarrange(unique_tax +labs(tag = 'A'), tpm_tax + labs(tag = 'B'), 
          common.legend = T, legend =  'bottom')
ggsave('out/ARGs/AMRf_unique_tpm_ARGs_taxonomy.png', dpi = 600)

## 
# Persistence of ARgs within an individual 
person_day <- amr_simplyfied %>%
  group_by(person) %>% 
  reframe(all_timepoints = n_distinct(day))

present <- amr_simplyfied %>% 
  group_by(person, ARG, day, Class) %>%
  reframe(PA = ifelse(sum(TPM) > 0, 1, 0)) %>% 
  group_by(person, ARG, Class) %>%
  reframe(timepoint_present = sum(PA), 
          timepoint_missing = sum(PA = 0)) %>% 
  left_join(person_day, by = 'person') %>%  
  mutate(prevalence = (timepoint_present/all_timepoints)*100) %>% 
  filter(prevalence > 0) %>%  
  group_by(person) %>% 
  mutate(all_args = n_distinct(ARG)) %>% 
  ungroup() %>% 
  group_by(person, Class, prevalence) %>% 
  reframe(n_args_class = n_distinct(ARG), 
          percent_args = (n_args_class/all_args) *100) 
  
present %>% 
  ggplot(aes(y = percent_args, x = prevalence, color = person)) +
  geom_point(size = 2) +
  geom_line(linewidth = 1.2) +
  #geom_smooth(se = F) +
  scale_x_continuous(breaks = seq(0, 100, by = 25)) +
  labs(y = '', x = '# timepoints an ARG was found', color = 'Individual') +
  facet_wrap(~Class, scales = 'free_x')
ggsave('out/ARGs/AMRf_args_present.png')


missing <- filter(present, prevalence < 90) 

amr %>%
  group_by(ARG, person, time_point, Class) %>%
  reframe(TPM = sum(TPM)) %>%
  filter(TPM > 0) %>%
  ggplot(aes(x = time_point, y = TPM, color = ARG)) +
  geom_line(show.legend = TRUE) +
  scale_y_log10() +
  facet_wrap(~person, scales = 'free_y')
ggsave('out/ARGs/tpm_time_ARGs.png')

amr %>%
  group_by(person, time_point, mechanism_resistence) %>%
  reframe(TPM = sum(TPM)) %>%
  filter(TPM > 0) %>%
  ggplot(aes(x = time_point, y = TPM, color = mechanism_resistence)) +
  geom_line(linewidth = 1) +
  scale_y_log10() +
  facet_wrap(~person, scales = 'free_y')
ggsave('out/ARGs/AMRf_persistence_mechanism_resistence.png')

# Correlation Shannons diveristy based on SGBs and TPM/unique ARGs
# At the level of Class of ARG
amr_simplyfied %>%
  group_by(person, day, Class) %>%
  reframe(TPM = sum(TPM)) %>%
  filter(TPM > 0) %>%
  ggplot(aes(x = day, y = TPM, color = Class)) +
  geom_line(linewidth =1.2) +
  scale_y_log10() +
  facet_wrap(~person, scales = 'free')+
  labs(x = 'Day', y = 'TPM [log10]', color = 'Class fo ARGs')
ggsave('out/ARGs/tpm_time_Class.png')

alpha <- readRDS('data/longitudinal_shotgun/alpha_diveristy.RDS')

tpm_amr <- amr_simplyfied %>%
  group_by(person, time_point, date, Class) %>%
  reframe(TPM = sum(TPM), 
          unique = n_distinct(ARG))

amr_alpha <- select(alpha, name, richness, shannon, person, date, time_point) %>% 
  right_join(tpm_amr, by = c('person', 'date', 'time_point'))

amr_alpha %>% 
  ggplot(aes(x = shannon, y = TPM)) +
  geom_smooth(se = F) +
  geom_point(size = 2, aes(color = Class)) +
  facet_wrap(~person, scales = 'free')
ggsave('out/ARGs/tpm_shannon.png')

