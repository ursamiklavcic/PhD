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

metadata <- read_csv2('data/metadata.csv') %>% filter(biota == 'bulk microbiota')

# Read in ARMfinder data 
amr_pre <- read_tsv('data/AMR/AMRfinder/amrfinder_results.tsv') %>%
  rename('contig_id' = `Contig id`, 'ARG' = `Element symbol`) %>%
  filter(`% Coverage of reference` > 80 & `% Identity to reference` > 50)

unique(amr_pre$ARG)

mechanisms <- read_csv2('data/AMR/AMRfinder/mechanism_resistence.csv')

# Run this lines only the first time
# Extract names of contigs for which ARGs were found
# orf_ids = amr_pre$contig_id
# writeLines(orf_ids, 'data/intermediate/AMRfinder_ORF.list')

# Before you run this script forward, make sure to downsize the 10.<>.mapcount file with Python script code/filter_AMRfinder_ORFs.sh
# bash filter_AMRfinder_ORFs.sh
contig_amr_coverage <- read.table('data/sqm_tables/19.filtered_AMRfinder.contigtable', header = TRUE, sep = '\t') %>%
  rename('contig_id' = 'Contig.ID') %>%
  select(-Bin.ID, -starts_with('Raw.'), -starts_with('TPM')) %>%
  pivot_longer(names_to = 'Group', values_to = 'coverage', cols = starts_with('Cov')) %>%
  #pivot_longer(names_to = 'Group2', values_to  = 'TPM', cols = starts_with('TPM')) %>%
  mutate(Group = str_sub(Group, 10, 14))

contig_amr_TPM <- read.table('data/sqm_tables/19.filtered_AMRfinder.contigtable', header = TRUE, sep = '\t') %>%
  rename('contig_id' = 'Contig.ID') %>%
  select(-Bin.ID, -starts_with('Raw.'), -starts_with('Cov')) %>%
  pivot_longer(names_to = 'Group', values_to  = 'TPM', cols = starts_with('TPM')) %>%
  mutate(Group = str_sub(Group, 5, 9)) 


amr <- left_join(amr_pre, contig_amr_coverage, by = 'contig_id', relationship = "many-to-many") %>%
  left_join(contig_amr_TPM, by = c('contig_id', 'Group', 'Disparity', 'GC.perc', 'Length', 'Num.genes', 'Tax')) %>%
  left_join(metadata, by ='Group', relationship = "many-to-many") %>%
  filter(!is.na(person)) %>%
  left_join(mechanisms, by = 'ARG')
   

amr %>% summarise(n_distinct(ARG))

amr %>% summarise(n_distinct(Class))

# unique classes 
unique <- amr %>%
  group_by(Class) %>%
  summarise(sum_unique = n_distinct(ARG)) %>%
  ggplot(aes(x = sum_unique, y = reorder(Class, sum_unique))) +
  geom_col() +
  labs(x = '# unique genes', y = 'Class of ARG')
unique
ggsave('out/ARGs/AMRf_class_unique.tiff')

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
ggsave('out/ARGs/AMRf_TPM_unique_class_ARGs_log.tiff')


ggarrange(unique, tpm, common.legend = TRUE)
ggsave('out/ARGs/AMRf_TPM_class_ARGs.tiff')

# mechanism of resistence 
amr %>%
  group_by(mechanism_resistence) %>%
  reframe(sum_TPM = sum(TPM),.groups = 'drop') %>%
  ggplot(aes(x = sum_TPM, y = reorder(mechanism_resistence, sum_TPM))) +
  geom_col() +
  scale_x_log10() +
  labs(y = 'Mechanism of resistence', x = 'log(TPM)') 
ggsave('out/ARGs/AMRf_mechanisms_tpm.tiff')

amr %>%
  group_by(mechanism_resistence) %>%
  reframe(unique = n_distinct(ARG),.groups = 'drop') %>%
  ggplot(aes(x = unique, y = reorder(mechanism_resistence, unique))) +
  geom_col() +
  scale_x_log10() +
  labs(y = 'Mechanism of resistence', x = '# unique ARGs') 
ggsave('out/ARGs/AMRf_mechanisms_unique.tiff')

# Number of different ARGs per individual through time 
amr %>% filter(TPM > 0) %>%
  group_by(person, time_point, Class) %>%
  summarise(sum_unique = n_distinct(ARG), .groups = 'drop') %>%
  filter(sum_unique > 0) %>%
  ggplot(aes(x = time_point, y = sum_unique, color  = Class)) +
  geom_point(size = 2) +
  geom_line(linewidth=1) +
  facet_wrap(~person, scales = 'free_y') +
  labs(y = '# unique genes', x = 'Time point', color = 'Class of ARG') 
ggsave('out/ARGs/AMRf_unique_class_ARG.tiff')

# TPM of different ARG classes through time 
amr %>% filter(TPM > 0 ) %>%
  group_by(person, time_point, Class) %>%
  summarise(sum_tpm = sum(TPM), .groups = 'drop') %>%
  filter(sum_tpm > 0) %>%
  ggplot(aes(x = time_point, y = sum_tpm, color  = Class)) +
  geom_point(size = 2) +
  geom_line(linewidth=1) +
  scale_y_log10() +
  facet_wrap(~person, scales = 'free_y') +
  labs(y = 'TPM', x = 'Time point', color = 'Class of ARG') 
ggsave('out/ARGs/AMRf_TPM_class_ARG_log.tiff')

# Taking into account events 
event_data <- amr %>%
  select(person, time_point, extremevent_type) %>%
  distinct() %>%
  filter(!is.na(extremevent_type)) %>% 
  mutate(xmin = time_point - 0.5, xmax = time_point + 0.5, ymin = -Inf,ymax = Inf)

amr %>%
  group_by(person, time_point, Class) %>%
  reframe(sum_TPM = sum(TPM),.groups = 'drop', extremevent_type) %>%
  ggplot(aes(x = time_point, y = sum_TPM)) +
  geom_rect(data = event_data, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = extremevent_type), inherit.aes = FALSE,
            alpha = 0.6) +
  scale_fill_manual(values = c('#d94343', '#d98e43', '#f1f011', '#0c9910', '#3472b7', 'white','#7934b7', '#b73485', '#0f5618')) +
  geom_point(size = 2, aes(color  = Class)) +
  geom_line(linewidth=1, aes(color  = Class)) +
  facet_wrap(~person, scales = 'free_y') +
  labs(y = 'log (TPM)', x = 'Time point', color = 'Class of ARG', fill = '') 
ggsave('out/ARGs/AMRf_extreme_event_TPM.tiff')

# What is the taxonomy of different ARG classes 
amr_tax <- amr %>%
  mutate(Phylum = str_extract(Tax, "p_[^;]+"),
         Phylum = substr(Phylum, 3, 25), 
         Kingdom = str_extract(Tax, 'k_[^;]+'), 
         Class2 = str_extract(Tax, 'c_[^;]+'))

amr_tax %>%
  ggplot(aes (x = Class, y = TPM, fill = Phylum)) +
  geom_boxplot() +
  scale_y_log10() +
  scale_fill_manual(values = c('#d94343', '#0c9910','#3472b7', '#b73485', '#f1f011', 'lightgrey' )) +
  theme_bw(base_size = 14) +
  facet_wrap(~Phylum, scales = 'free_y', ncol = 2) +
  theme(legend.position = 'bottom') +
  coord_flip()
ggsave('out/ARGs/AMRf_taxonomy_class_TPM.tiff')


# Hypothesis 3: Does an increase of ARGs correlate with decrease in diversity

##
## Alpha diversity 
##

# Is this in any way correlated with alpha diversity? 

otutabEM <- readRDS('~/projects/longitudinal_amplicons/data/r_data/otutabEM.RDS')
richnessEM = estimateR(otutabEM) # observed richness and Chao1
evennessEM = diversity(otutabEM)/log(specnumber(otutabEM)) # evenness index
shannonEM = diversity(otutabEM, index = 'shannon')

# Join all calculations and metadata
alpha_meta = as_tibble(as.list(evennessEM)) %>% pivot_longer(names_to = 'Group', values_to = 'evenness', cols = starts_with(c('M', 'S'))) %>%
  left_join(t(richnessEM) %>% as.data.frame() %>% rownames_to_column('Group'), by='Group') %>%
  left_join(as_tibble(as.list(shannonEM)) %>% pivot_longer(names_to = 'Group', values_to = 'shannon', cols = starts_with(c('M', 'S')))) %>%
  left_join(metadata, by='Group') %>%
  mutate(person2 = person) 

# Function to calculate correlation values for each ARGs class and the shannon diversity
amr_alpha <- amr %>%
  group_by(Group, person, time_point, Class) %>%
  filter(TPM > 0) %>%
  reframe(sum_tpm = sum(TPM), 
          sum_unique = n_distinct(ARG)) %>%
  left_join(alpha_meta, by = c('Group', 'person', 'time_point'))

amr_alpha %>% 
  ggplot(aes(x = shannon, y = log10(sum_tpm))) +
  geom_point(mapping = aes(color = person), size = 2) +
  geom_smooth(method = 'lm') +
  stat_cor() +
  facet_wrap(~Class, scales = 'free') +
  labs(x = 'Shannon diveristy index', y = 'log (TPM)', color = 'Individual') +
  theme(legend.position = 'bottom') +
  guides(color = guide_legend(nrow = 1))
ggsave('out/ARGs/AMRf_corr_tpm.tiff')


amr_alpha %>% 
  ggplot(aes(x = shannon, y = sum_unique)) +
  geom_point(mapping = aes(color = person), size = 2) +
  geom_smooth(method = 'lm') +
  stat_cor() +
  scale_y_continuous(label = scales::comma) +
  facet_wrap(~Class, scales = 'free_y') +
  labs(x = 'Shannon diveristy index', y = '# ARGs', color = 'Individual') +
  theme(legend.position = 'bottom') +
  guides(color = guide_legend(nrow = 1))
ggsave('out/ARGs/AMRf_corr_unique.tiff')

# Correlations calculation 
results = data.frame()
for (j in unique(amr$Class)) {
  x = amr %>% 
    filter(Class == j & TPM > 0 & !is.na(TPM)) %>%
    group_by(Group, person, time_point, Class) %>%
    reframe(sum_tpm = sum(TPM), 
            sum_unique = n_distinct(ARG)) %>%
    left_join(alpha_meta, by = c('Group', 'person', 'time_point'))
  
  corr_tpm <- cor.test(x$sum_tpm, x$shannon, method = 'pearson')
  corr_unique <- cor.test(x$sum_unique, x$shannon, method = 'pearson')
  
  results = rbind(results, data.frame(
    Class = j, 
    corr_tpm = corr_tpm$estimate, 
    p_tpm = corr_tpm$p.value, 
    corr_uniq = corr_unique$estimate, 
    p_uniq = corr_unique$p.value))
  
}

results

## Correlations between unique number of ARGs and shannon 
filter(results, !is.na(corr_uniq)) %>% 
  mutate(biota = 'Microbiota',
         signif_label = case_when(p_uniq <= 0.001 ~ "***", p_uniq <= 0.01 ~ "**",
           p_uniq <= 0.05 ~ "*", TRUE ~ ""), 
         label = paste0(sprintf("%.3f", corr_uniq), ' ', signif_label)) %>%
  ggplot(aes(x = biota, y = Class, fill = corr_uniq)) +
  geom_tile() +
  geom_text(aes(label = label), color = 'black', size = 4) +
  scale_fill_gradient2(low = "#3472b7", mid = "white", high = "#0c9910", midpoint = 0) +
  labs(
    caption = "Correlation between the number of unique ARGs \n and Shannon's diversity index based on OTUs",
    y = '', x = '', fill = "Correlation Coefficient"
  ) +
  theme_bw(base_size = 12)

ggsave('out/ARGs/AMFf_corr_unique.tiff')

# between TPM to ARGs and shannon
filter(results, !is.na(corr_uniq)) %>% 
  mutate(biota = 'Microbiota',
         signif_label = case_when(p_tpm <= 0.001 ~ "***", p_tpm <= 0.01 ~ "**",
                                  p_tpm <= 0.05 ~ "*", TRUE ~ ""), 
         label = paste0(sprintf("%.3f", corr_tpm), ' ', signif_label)) %>%
  ggplot(aes(x = biota, y = Class, fill = corr_tpm)) +
  geom_tile() +
  geom_text(aes(label = label), color = 'black', size = 4) +
  scale_fill_gradient2(low = "#3472b7", mid = "white", high = "#0c9910", midpoint = 0) +
  labs(
    caption = "Correlation between the TPM of ARGs \n and Shannon's diversity index based on OTUs",
    y = '', x = '', fill = "Correlation Coefficient"
  ) +
  theme_bw(base_size = 12)

ggsave('out/ARGs/AMFf_corr_tpm.tiff')

## 
# Persistence of ARgs within an individual 
present <- amr %>%
  filter(time_point < 13) %>%
  group_by(person, time_point, ARG, Class) %>%
  reframe(sumTPM = sum(TPM), 
          PA = ifelse(sumTPM > 0, 1, 0)) %>%
  group_by(person, ARG, Class) %>%
  reframe(all_timepoints = n(), 
          timepoints_present = sum(PA == 1), 
          timepoints_missing = sum(PA == 0))
  
present %>% 
  ggplot(aes( x= Class, y = timepoints_present, fill = Class)) +
  geom_violin(draw_quantiles = c( 0.5)) +
  scale_y_continuous(breaks = c(1:12)) +
  labs(x = '', y = '# timepoints an ARG was found', fill = '') +
  theme(legend.position = 'none') +
  coord_flip()
ggsave('out/ARGs(AMRf_args_present.tiff')


missing <- filter(present, timepoints_present < 11) 

amr %>%
  group_by(ARG, person, time_point, Class) %>%
  reframe(TPM = sum(TPM)) %>%
  filter(TPM > 0) %>%
  ggplot(aes(x = time_point, y = TPM, color = ARG)) +
  geom_line(show.legend = TRUE) +
  scale_y_log10() +
  facet_wrap(~person, scales = 'free_y')
ggsave('out/ARGs/tpm_time_ARGs.tiff')

amr %>%
  group_by(person, time_point, mechanism_resistence) %>%
  reframe(TPM = sum(TPM)) %>%
  filter(TPM > 0) %>%
  ggplot(aes(x = time_point, y = TPM, color = mechanism_resistence)) +
  geom_line(linewidth = 1) +
  scale_y_log10() +
  facet_wrap(~person, scales = 'free_y')

# At the level of Class of ARG
amr %>%
  group_by(person, time_point, Class) %>%
  reframe(TPM = sum(TPM)) %>%
  filter(TPM > 0) %>%
  ggplot(aes(x = time_point, y = TPM, color = Class)) +
  geom_line() +
  scale_y_log10() +
  facet_wrap(~person, scales = 'free_y')
ggsave('out/ARGs/tpm_time_Class.tiff')

# Prevalence of ARGs, mechanisms, classes of ARGs 

# Third plot % OTUs on y, x 0 prevalence % 
prevalence <- amr %>%
  mutate(PA = ifelse(TPM > 0, 1, 0)) %>%
  group_by(person, ARG, Class, mechanism_resistence) %>%
  reframe(all_timepoints = n(), 
          timepoints_present = sum(PA == 1), 
          timepoints_missing = sum(PA == 0)) %>%
  # OTU had to be present in at least 50% of all samples from 1 individual! 
  # Remove 'singletons'
  #filter(timepoints_present > 1) %>%
  mutate(prevalence = (timepoints_present/all_timepoints)*100) 

prevalence_ARGs <- prevalence %>%
  group_by(person) %>%
  mutate(all_args = n_distinct(ARG)) %>%
  ungroup() %>%
  group_by(person, prevalence) %>%
  reframe(no_args = n_distinct(ARG), 
          per_args = (no_args/all_args) *100) 

prevalence_ARGs %>% 
  ggplot(aes(x = prevalence, y = per_args, color = person)) +
  geom_line(linewidth = 1)


prevalence_class <-prevalence %>%
  group_by(person, Class) %>%
  mutate(all_args = n_distinct(ARG)) %>%
  ungroup() %>%
  group_by(person, Class, prevalence) %>%
  reframe(no_args = n_distinct(ARG), 
          per_args = (no_args/all_args) *100) 

prevalence_class %>%
  ggplot(aes(x = prevalence, y = per_args, color = Class)) +
  geom_line(linewidth = 1) +
  facet_grid(~person)
