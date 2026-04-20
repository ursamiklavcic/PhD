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
theme_set(theme_bw(base_size = 12) +
            theme(plot.title   = element_text(size = 12),
                  axis.title   = element_text(size = 12),
                  axis.text    = element_text(size = 12)))

metadata <- read_csv2('~/projects/longitudinal_shotgun/data/metadata.csv') 

# Read in ARMfinder data 
amr_pre <- read_tsv('~/projects/longitudinal_shotgun/data/AMR/AMRfinder/amrfinder_results.tsv') %>%
  rename('contig_id' = `Contig id`, 'ARG' = `Element symbol`) %>%
  filter(`% Coverage of reference` > 80 & `% Identity to reference` > 50) %>%  
  filter(Type == 'AMR') %>% 
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
          nrow = 2)
ggsave('out/ARGs/unique_tpm_args_mechanism.svg', dpi = 600)


# TPM of different ARG classes through time 
# Taking into account events 
event_data <- read.table('data/extreme_event_data.csv', sep = ',', header = TRUE) 

tpm_time <- amr_simplyfied %>%
  group_by(person, day, Class) %>%
  reframe(sum_TPM = sum(TPM), extremevent_type) %>%
  ggplot(aes(x = day, y = sum_TPM)) +
  geom_rect(data = event_data, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = extremevent_type), inherit.aes = FALSE,
            alpha = 0.6) +
  scale_fill_manual(values = c('white','#d94343', '#d98e43', '#f1f011', '#0c9910', '#3472b7', '#7934b7', '#b73485', '#0f5618')) +
  geom_point(size = 2, aes(color  = Class)) +
  geom_line(linewidth = 1, aes(color  = Class)) +
  facet_wrap(~person, scales = 'free_y') +
  labs(y = 'ARG abundance [log10(TPM)]', x = 'Day', color = 'Class of ARG', fill = 'Event') 
tpm_time
ggsave('out/ARGs/AMRf_extreme_event_TPM.png')

# All ARGs together 
amr_simplyfied %>%
  group_by(person, day, ) %>%
  reframe(sum_TPM = sum(TPM), extremevent_type) %>%
  ggplot(aes(x = day, y = sum_TPM)) +
  geom_rect(data = event_data, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = extremevent_type), inherit.aes = FALSE,
            alpha = 0.6) +
  scale_fill_manual(values = c('white','#d94343', '#d98e43', '#f1f011', '#0c9910', '#3472b7', '#7934b7', '#b73485', '#0f5618')) +
  geom_point(size = 2) +
  geom_line(linewidth = 1) +
  facet_wrap(~person) +
  labs(y = 'ARG abundance [log10(TPM)]', x = 'Day', fill = 'Event') 
ggsave('out/ARGs/all_ARGs_time.png', dpi=600)

# Number of different ARGs per individual through time 
unique_time <- amr_simplyfied %>% filter(TPM > 0) %>%
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
  labs(y = 'ARG diversity\n[number of unique genes]', x = 'Day', color = 'Class of ARG', fill = 'Event') 
unique_time
ggsave('out/ARGs/AMRf_unique_class_ARG.png')

ggarrange(unique_time + labs(tag = 'A'), 
          tpm_time + labs(tag = 'B'), 
          common.legend = T, 
          legend = 'bottom', 
          ncol = 1)
ggsave('out/ARGs/unique_tpm_time_up.svg', dpi=600)

# Unique all 
amr_simplyfied %>% filter(TPM > 0) %>%
  group_by(person, day) %>%
  summarise(sum_unique = n_distinct(ARG), .groups = 'drop') %>%
  filter(sum_unique > 0) %>%
  ggplot(aes(x = day, y = sum_unique)) +
  geom_rect(data = event_data, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = extremevent_type), inherit.aes = FALSE,
            alpha = 0.6) +
  scale_fill_manual(values = c('white','#d94343', '#d98e43', '#f1f011', '#0c9910', '#3472b7', '#7934b7', '#b73485', '#0f5618')) +
  geom_point(size = 2) +
  geom_line(linewidth=1) +
  facet_wrap(~person, scales = 'free_y') +
  labs(y = '# unique genes', x = 'Day', fill = 'Event type')
ggsave('out/ARGs/unique_ARGs_all.png', dpi=600)

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
unique_tax
ggsave('out/ARGs/AMRf_taxonomy_class_unique_simple.png')

# Which phylum has the highest abundance of ARGs? 
tpm_tax <- amr_tax %>%
  group_by(Class, Phylum) %>%  
  reframe(TPM = sum(TPM)) %>%  
  ggplot(aes (x = TPM, y = Phylum, fill = Class)) +
  geom_col() +
  scale_x_log10() +
  #scale_fill_manual(values = c('#d94343', '#0c9910','#3472b7', '#b73485', '#f1f011', 'lightgrey' )) +
  theme(legend.position = 'bottom', axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  labs(y = '', x = 'TPM [log10]')
tpm_tax
ggsave('out/ARGs/AMRf_taxonomy_class_TPM_simple_unique.png')

ggarrange(unique_tax +labs(tag = 'A'), tpm_tax + labs(tag = 'B'), 
          common.legend = T, legend =  'bottom')
ggsave('out/ARGs/AMRf_unique_tpm_ARGs_taxonomy.png', dpi = 600)

# Mechanisms through time
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
alpha <- readRDS('data/longitudinal_shotgun/alpha_diveristy.RDS') %>%  
  filter(biota == 'bulk microbiota')

tpm_unique_amr <- amr_simplyfied %>%
  filter(TPM > 0) %>% 
  group_by(Group, person, Class) %>%
  reframe(TPM = sum(TPM), 
          unique = n_distinct(ARG))

amr_alpha <- select(alpha, name, richness, shannon, evenness, date, day) %>% 
  full_join(tpm_unique_amr, by = join_by('name' == 'Group'))

# Correlation linear-mixed effect model for correlation between Shannon diveristy and ARG diveristy and/or abundance 
amr_alpha

library(lmerTest)

lme_simple <-lmerTest::lmer(data = amr_alpha, formula = unique ~ shannon + day + (1 | person))
summary(lme_simple)

# Does this positive interaction change with time? 
lme_time <- lmerTest::lmer(data = amr_alpha, formula = unique ~ shannon * day + (1 | person)) # NO time does not have an effect! 
summary(lme_time)

# Per class? 
class1_lme <- lmerTest::lmer(data = amr_alpha, formula = unique ~ shannon + (1 | person) + (1 + shannon | Class))
summary(class1_lme)

# What if Class is fixed effect ? 
class_lme <- lmerTest::lmer(data = amr_alpha, formula = unique ~ shannon * Class + (1 | person))
summary(class_lme)


# Plot 
coefs <- summary(class_lme)$coefficients
labels <- coefs[grep("shannon:Class", rownames(coefs)), "Pr(>|t|)"] %>% 
  as.data.frame() %>% 
  rownames_to_column('Class') %>% 
  mutate(Class = substr(Class, 14, 100)) %>% 
  rbind(data.frame(Class = 'AMINOGLYCOSIDE', . = 0.4414))

ggplot(amr_alpha, aes(x = shannon, y = unique, color = person)) +
  geom_point(alpha = 0.7, size = 2) +
  geom_smooth(method = "lm", se = TRUE, color = "black") +
  geom_text(data = labels, aes(x = 4, y = 2.5, label = paste0('p = ',signif(., 2))), inherit.aes = F) +
  facet_wrap(~ Class, scales = "free_y") +
  labs(x = "Shannon's diversity index",
    y = "ARG diversity [# unique genes]", color = 'Individual')
ggsave('out/ARGs/lme_corr_simple.svg', dpi=400)

# testing abundance ARGs + bacterial diveristy! 
classTPM_lme <- lmerTest::lmer(data = amr_alpha, formula = log10(TPM) ~ shannon * Class + (1 | person))
summary(classTPM_lme, correlation=TRUE)


# # Matrix of correlations between Shannon and abundance of ARGs 
# cor_data <- amr_alpha %>%
#   group_by(person, Class) %>%
#   reframe(cor_shannon_tpm = cor.test(shannon, TPM, method = "pearson")$p.value,
#           cor_shannon_unique = cor.test(shannon, unique, method = "pearson")$p.value)
# 
# ggplot(amr_alpha, aes(x = shannon, y = TPM)) +
#   geom_point(aes(color = person)) +
#   geom_smooth(method = 'lm') +
#   stat_cor(label.y.npc = 'top', method = 'spearman') +
#   facet_wrap(~Class, scales = 'free_y') +
#   labs(x = 'Shannon', y = 'TPM [log10]', color = 'Individual')
# ggsave('out/ARGs/correlation_abundanceARG_shannon.svg', dpi = 600)
# 
# ggplot(cor_data, aes(x = person, y = Class, fill = cor_shannon_tpm)) +
#   geom_tile() +
#   geom_text(aes(label = sprintf("%.3f", cor_shannon_tpm))) +
#   scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0.1, limits = c(0, 1)) + 
#   labs(fill = "Pearson\n correlation", x = "Individual", y = '', 
#        title = 'Correlation between Shannon diversity and abundance of an ARG class') 
# ggsave('out/ARGs/correlation_matrix_abundanceARG_shannon.png', dpi=600)
# 
# # cor.test(filter(amr_alpha, person == 'A')$shannon, filter(amr_alpha, person == 'A')$TPM, method = 'pearson')$p.value
# 
# # Matrix shannon and diveristy of ARGs 
# ggplot(amr_alpha, aes(x = shannon, y = unique)) +
#   geom_point(aes(color = person)) +
#   geom_smooth(method = 'lm') +
#   stat_cor(label.y.npc = 'bottom', method = 'spearman') +
#   scale_y_log10() +
#   facet_wrap(~Class) +
#   labs(x = 'Shannon', y = '# unique ARGs', color = 'Individual')
# ggsave('out/ARGs/correlation_diversityARG_shannon.svg', dpi=600)
# 
# ggplot(cor_data, aes(x = person, y = Class, fill = cor_shannon_unique)) +
#   geom_tile() +
#   geom_text(aes(label = sprintf("%.3f", cor_shannon_unique))) +
#   scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0.1, limits = c(0, 1)) + 
#   labs(fill = "Pearson\n p-value", x = "Individual", y = '', 
#        title = 'Correlation between Shannon diversity and ARG diveristy') 
# ggsave('out/ARGs/correlation_matrix_diversityARG_shannon.png', dpi=600)
# 
# # For all ARGs together = inflating statistics! 
# amr_group <- amr_simplyfied %>%
#   filter(TPM > 0) %>% 
#   group_by(Group, person) %>%
#   reframe(TPM = sum(TPM), 
#           unique = n_distinct(ARG))
# 
# alpha_amr_2 <- select(alpha, name, richness, shannon, evenness, date) %>% 
#   full_join(amr_group, by = join_by('name' == 'Group'))
# 
# ggplot(alpha_amr_2, aes(x = shannon, y = TPM)) +
#   geom_point(aes(color = person)) +
#   geom_smooth(method = 'lm') +
#   stat_cor(label.y.npc = 'top', method = 'pearson') +
#   facet_wrap(~person, scales = 'free') +
#   labs(x = 'Shannon', y = 'TPM of all ARGs [log10]', color = 'Individual')
# ggsave('out/ARGs/all_ARG_abundance_shannon.png', dpi = 600)
# 
# ggplot(alpha_amr_2, aes(x = shannon, y = unique)) +
#   geom_point(aes(color = person)) +
#   geom_smooth(method = 'lm') +
#   stat_cor(label.y.npc = 'top', method = 'pearson') +
#   facet_wrap(~person, scales = 'free') +
#   labs(x = 'Shannon', y = '# unique ARGs', color = 'Individual')
# ggsave('out/ARGs/all_ARG_diversity_shannon.png', dpi = 600)

# CTX-M1
amr %>% filter(ARG == 'blaCTX-M-1') %>% 
  group_by(day, person, ARG) %>% 
  reframe(TPM = sum(TPM)) %>% 
  filter(TPM > 0.5) %>% 
  ggplot(aes(x = day, y = TPM, color = person)) +
  geom_point(size = 3) +
  scale_y_log10() +
  theme(legend.position = 'none') +
  labs(x = 'Day', y = 'TPM [log10]', title = 'Individual I')
ggsave('out/ARGs/ctx-m1.png')  

unique(filter(amr, ARG == 'blaCTX-M-1', TPM > 0.5)$Tax)

# Van resistence 
van <- amr %>%  
  mutate(van = case_when(str_detect(ARG, 'van') ~ 'van', 
                         TRUE ~ 'not van')) %>% 
  filter(van == 'van') 

unique(filter(van, TPM > 0.5)$Tax)

van %>% 
  group_by(day, person, ARG, Tax) %>% 
  reframe(TPM = sum(TPM)) %>% 
  filter(TPM > 0.5) %>%  
  ggplot(aes(x = day, y = TPM, color = Tax)) +
  geom_point(size = 2) +
  geom_line(linewidth=1.2) +
  scale_y_log10() +
  facet_wrap(ARG~person) +
  labs(x = 'Day')
ggsave('out/ARGs/van.png')


## 
# Persistence of ARgs within an individual 
# person_day <- amr %>%
#   group_by(person) %>% 
#   reframe(all_timepoints = n_distinct(day))
# 
# present <- amr %>% 
#   group_by(person, ARG, day, Class) %>%
#   reframe(PA = ifelse(sum(TPM) > 0, 1, 0)) %>% 
#   group_by(person, ARG, Class) %>%
#   reframe(timepoint_present = sum(PA), 
#           timepoint_missing = sum(PA = 0)) %>% 
#   left_join(person_day, by = 'person') %>%  
#   mutate(prevalence = (timepoint_present/all_timepoints)*100) %>% 
#   filter(prevalence > 0) %>%  
#   group_by(person) %>% 
#   mutate(all_args = n_distinct(ARG)) %>% 
#   ungroup() %>% 
#   group_by(person, Class, prevalence) %>% 
#   reframe(n_args_class = n_distinct(ARG), 
#           percent_args = (n_args_class/all_args) *100) 
# 
# present %>% 
#   filter(!Class %in% c('FOSFOMYCIN', 'LINCOSAMIDE/MACROLIDE', 'LINCOSAMIDE/STREPTOGRAMIN', 
#                        'MACROLIDE/STREPTOGRAMIN', 'PHENICOL/LINCOSAMIDE/OXAZOLIDINONE/PLEUROMUTILIN/STREPTOGRAMIN', 
#                        'QUINOLONE', 'SULFONAMIDE', 'PHENICOL'), n_args_class > 1) %>% 
#   ggplot(aes(y = n_args_class, x = prevalence, color = Class)) +
#   geom_point(size = 3) +
#   geom_line(linewidth = 2, alpha = .7) +
#   labs(y = '', x = '% timepoints an ARG was found', color = 'Individual') +
#   facet_wrap(~person, scales = 'free')
# ggsave('out/ARGs/AMRf_args_present.png')
# 
# 
# present %>% filter(n_args_class > 1)

# Heatmap genes x time for each person
# Only ARGs present in more than 1 timepoint and not always
amr %>%  
  group_by(person, ARG, day, Class) %>%
  reframe(PA = ifelse(sum(TPM) > 0, 1, 0)) %>% 
  group_by(person, ARG, Class)  %>%  
  mutate(sum = sum(PA)) %>% 
  filter(sum > 1) %>% 
  filter(sum < 11) %>% 
  ggplot(aes(x = as.factor(day), y = ARG, fill = as.factor(PA))) +
  geom_tile() +
  scale_fill_manual(values = c('white', 'blue')) +
  facet_wrap(~person, scales = 'free') +
  theme(axis.text.y = element_text(size = 6))
ggsave('out/ARGs/more_than_once_but_not_always.png')


# Always present ARGs 
amr %>%  
  group_by(person, ARG, day, Class) %>%
  reframe(PA = ifelse(sum(TPM) > 0, 1, 0)) %>% 
  group_by(person, ARG, Class)  %>%  
  mutate(sum = sum(PA)) %>% 
  filter(sum > 11) %>% 
  ggplot(aes(x = as.factor(day), y = ARG, fill = as.factor(PA))) +
  geom_tile() +
  scale_fill_manual(values = c('white', 'blue')) +
  facet_wrap(~person, scales = 'free') +
  theme(axis.text.y = element_text(size = 6))
ggsave('out/ARGs/always_present.png')

# ARGs persisitence
amr_persistence <- amr %>%  
  group_by(person, ARG, day, Class, TPM) %>%
  reframe(PA = ifelse(sum(TPM) > 0, 1, 0)) %>% 
  group_by(person, ARG, Class)  %>%  
  mutate(sum = sum(PA)) 

# Diiferent ARGS per person
amr %>% filter(TPM > 0) %>% 
  group_by(person) %>%  
  reframe(len = n_distinct(ARG))

# How many ARGs are always present ? 
amr_persistence %>% 
  filter(sum > 11) %>% 
  group_by(person) %>%  
  reframe(len = n_distinct(ARG))
# Only once
amr_persistence %>% 
  filter(sum < 2) %>% 
  group_by(person) %>%  
  reframe(len = n_distinct(ARG))

# Table for thesis 
table_amrs <- amr %>% 
  group_by(Class, ARG, person, day) %>% 
  reframe(PA = ifelse(TPM > 0, 1, 0)) %>%  
  group_by(Class, ARG, person) %>%  
  reframe(sum = sum(PA)) %>% 
  filter(sum > 10) %>%  
  group_by(Class, ARG) %>% 
  reframe(n = n_distinct(person))
write_csv2(table_amrs, 'out/ARGs/table.csv')
