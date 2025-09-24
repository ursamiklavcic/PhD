# Analysis of V3V4 16S rRNA data for comparison of ethanol/temperature shock and the usage of EMA
# for DNA isolations of sporeforming bacteria in human stool samples
# OTUs were constructed in mothur, code availabile HOPC/bin/mothur.script 

library(dplyr)
library(tidyr)
library(tibble)
library(vegan)
library(stringr)
library(readr)
library(ggplot2)
library(scales)
library(glue)

set.seed(96)
theme_set(theme_bw(base_size = 12))

# Exploration
shared = read.table('etoh_t_comparison/data/mothur/stability.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.opti_mcc.shared', 
                    sep = '\t', header = TRUE) %>%
  select(Group, starts_with('Otu')) %>%
  pivot_longer(-Group)

# Distribution of reads per sample 
shared %>%
  group_by(Group) %>%
  summarize(ReadsPerSample=sum(value)) %>%
  ggplot(aes(x=ReadsPerSample)) +
  geom_histogram() +
  scale_x_continuous(breaks = seq(0, 80000, by=10000)) +
  #scale_y_continuous(breaks = seq(0,40, by=5)) +
  labs(x = 'Reads per sample', y = 'Number of samples')
ggsave('etoh_t_comparison/plots/reads_per_sample.png', dpi=600)

# How min/max/mean/sum reads per sample/all samples
info1 = shared %>%
  group_by(Group) %>%
  summarize(ReadsPerSample=sum(value))
min(info1$ReadsPerSample) # 3
max(info1$ReadsPerSample) # 85784
mean(info1$ReadsPerSample) # 25006.03
median(info1$ReadsPerSample) # 23534
sum(info1$ReadsPerSample) # 4000965

# Distribution of OTUs
shared %>%
  group_by(name) %>%
  summarize(OTUabundance=sum(value)) %>%
  ggplot(aes(x=OTUabundance)) +
  geom_histogram(breaks=seq(0, 75, by =1)) +
  labs(x = 'Abundance of OTU', y= 'Number of OTUs')
ggsave('etoh_t_comparison/plots/reads_per_OTU.png', dpi=600)

# How min/max/mean reads per OTU
info2 = shared %>%
  group_by(name) %>%
  summarize(OTUabundance=sum(value))
min(info2$OTUabundance) # 1
max(info2$OTUabundance) # 550584
mean(info2$OTUabundance) # 249.3435
median(info2$OTUabundance) # 1
sum(info2$OTUabundance) # 4000965 = Total number of redas in the dataset before removing anything! 

reads_per_OTU = shared %>%
  group_by(name) %>%
  summarize(OTUabundance=sum(value))
# Number of all OTUs 
length(unique(shared$name)) # 16046

# How many OTUs have less thaN 10 reads
sum(reads_per_OTU$OTUabundance < 10) # 15100
# How many reads do they contain
reads_per_OTU %>%
  filter(OTUabundance < 10) %>%
  summarise(sum(OTUabundance)) # 22 486

# what percentage is this from all reads ?
sum(reads_per_OTU$OTUabundance < 10)/sum(shared$value) *100 #  0.4%

# Min number of reads in sample and min sequences per OTU as determined in the exploration analysis
reads_per_sample = 10000
min_seqs_per_otu = 10

shared_pre = shared %>%
  group_by(Group) %>%
  # Count the number of reads in each sample
  mutate(sum_sample = sum(value)) %>%
  # and remove all from microbiota and sporobiota, that have less than 100 000 
  filter(sum_sample > reads_per_sample) %>%
  ungroup() %>%
  group_by(name) %>%
  # Count the number of reads each OTU has
  mutate(sum_otus = sum(value)) %>%
  # Remove OTUs that have less than 0,000001% reads in total
  filter(sum_otus > min_seqs_per_otu) %>%
  ungroup() %>%
  select(-sum_otus, -sum_sample)

# 
otutab_pre = shared_pre %>%
  pivot_wider(names_from = 'name', values_from = 'value', values_fill = 0) %>%
  column_to_rownames('Group')

# Rarefy the data once - as we sequenced so deep that for the first analysis this is not crucial !
otutab = rrarefy(otutab_pre, sample=reads_per_sample)
saveRDS(otutab, 'etoh_t_comparison/data/r_data/otutab.RDS')

# Extract OTUs that are present rarefied table 
otu_names = as.data.frame(otutab) %>% colnames() 

# Import taxonomy table
taxtab = read_tsv('etoh_t_comparison/data/mothur/stability.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.opti_mcc.0.03.cons.taxonomy') %>%
  filter(OTU %in% otu_names) %>%
  select(name = "OTU", taxonomy = "Taxonomy") %>%
  mutate(taxonomy = str_replace_all(taxonomy, "\\\\|\\\"|\\(\\d+\\)", ""),
         taxonomy = str_replace(taxonomy, ";$", "")) %>%
  separate(taxonomy, into=c("Domain", "Phylum", "Class", "Order", "Family", "Genus"),
           sep=";") %>%
  mutate(Phylum = ifelse(Phylum %in% c('Firmicutes', 'Bacteroidetes', 'Actinobacteria', 'Proteobacteria',
                                       'Verrucomicrobia','Bacteria_unclassified'), Phylum, 'Other')) %>%
  mutate(Phylum = case_when(
    Phylum == 'Firmicutes' ~ 'Bacillota',
    Phylum == 'Bacteroidetes' ~ 'Bacteroidota',
    Phylum == 'Actinobacteria' ~ 'Actinomycetota',
    Phylum == 'Proteobacteria' ~ 'Pseudomonadota',
    Phylum == 'Bacteria_unclassified' ~ 'unclassified Bacteria',
    Phylum == 'Verrucomicrobia' ~ 'Verrucomicrobiota', 
    #Phylum == 'Lentisphaerae' ~ 'Lentisphaerota', 
    #Phylum == 'Fusobacteria' ~ 'Fusobacteriota', 
    # Phylum == 'Synergistetes' ~ 'Synergistota',
    TRUE ~ Phylum ))
saveRDS(taxtab, 'etoh_t_comparison/data/r_data/taxtab.RDS') 
  
unique(taxtab$Phylum)

# Import metadata
metadata <- as_tibble(read.csv('etoh_t_comparison/data/metadata.csv', sep=',')) %>%
  filter(Group %in% rownames(otutab))

# remove unnecessary 
rm(info1)
rm(info2)
rm(otutab_pre)
rm(min_seqs_per_otu)
rm(otu_names)
rm(reads_per_sample)
rm(shared_pre)
rm(reads_per_OTU)
rm(shared)

# Plots

col_phylum = c('#1F77B4', '#FF7F0E',  '#2CA02C',  '#D62728', '#9467BD', '#8C564B', '#f4d03f')
col_shock = c('#a569bd', '#27ae60', '#28a9d8')

################# DNA CONCENTRATIONS ########################
conc = metadata %>% 
  group_by(treatment) %>%
  mutate(meanConc = mean(DNAconcentration)) 

conc %>% ggplot(aes(x=DNAconcentration, y=treatment, fill = treatment)) +
  geom_boxplot() +
  geom_point(color = 'black') +
  scale_y_discrete(labels=c('Ethanol shock (100%) + cultivation', 'Ethanol shock (70%) + cultivation','Heat shock (30 min) + cultivation', 'Heat shock (60 min) + cultivation',
                            'Ethanol shock + EMA 10 min', 'Heat shock + EMA 10 min', 'Ethanol shock + EMA 15 min', 'Heat shock + EMA 15 min', 'Ethanol shock + EMA 5 min',
                            'Heat shock + EMA 5 min', 'Microbiota', 'Ethanol shock + wash',  'Heat shock + wash')) +
  ylab('Treatment') +
  xlab('DNA concentration (ng/μl)') +
  theme(legend.position="none")
ggsave('etoh_t_comparison/plots/DNA_concentration_treatment.png', dpi=600)


# CFU counts 
cfu <- as_tibble(read.csv('etoh_t_comparison/data/CFU.csv', sep=',')) %>%
  mutate(bile_acids = ifelse(bile_acids == 'yes', 'With bile acids', 'Without bile acids'), 
         cultivation_media = ifelse(cultivation_media == 'liquid', 'Liquid media', 'Solid media'), 
         treatment = ifelse(treatment == 'cult_70', 'Ethanol shock (70%)', 
                            ifelse(treatment == 'cult_100', 'Ethanol shock (100%)', 
                                   ifelse(treatment == 'cult_T30', 'Heat shock (30 min)', 'Heat shock (60 min)'))))

cfu %>%
  group_by(shock, treatment, cultivation_day, cultivation_media, bile_acids) %>%
  reframe(mean_CFU = mean(CFU, na.rm = TRUE)) %>%
  ungroup() %>%  
  ggplot(aes(x = cultivation_day, y=mean_CFU, color=treatment)) +
  geom_point(size = 3, position = position_dodge(width=0.2)) +
  geom_line(aes(linetype = bile_acids), linewidth = 2, position = position_dodge(width=0.2)) +
  facet_grid(~cultivation_media) +
  labs(x = 'Cultivation day', y = 'CFU', color = 'Treatment', linetype = '') +
  theme(legend.position = 'right')
ggsave("etoh_t_comparison/plots/CFU.png", dpi=600)

# Prepare long format of OTUtable 
otu_long <- as.data.frame(otutab) %>%
  rownames_to_column('Group') %>%
  left_join(metadata, by = 'Group') %>%
  pivot_longer(names_to = 'name', values_to = 'value', cols = starts_with('Otu')) %>%
  group_by(Group) %>%
  mutate(rel_abund = value/sum(value), 
         PA = ifelse(value > 0, 1, 0)) %>%
  ungroup() %>%
  left_join(taxtab, by = 'name')

otu_long$Phylum <- factor(otu_long$Phylum, levels = c('Bacillota', 'Bacteroidota', 'Actinomycetota', 
                                                      'Pseudomonadota','Verrucomicrobiota', 'unclassified Bacteria', 'Other'))
###
## Cultivation experiment heat VS ethanol shock
###

otu_cultivation <- filter(otu_long, cultivation %in% c('cultivationSolid', 'cultivationLiquid') | Group %in% c('UM', 'TZ')) %>%
  filter(!is.na(Phylum)) %>%  
  mutate(bile_acids = ifelse(bile_acids == 'yes', 'With bile acids', 
                             ifelse(bile_acids == 'no', 'Without bile acids', 'Stool sample')), 
         cultivation_media = ifelse(cultivation_media == 'liquid', 'Liquid media',
                                    ifelse(cultivation_media == 'solid', 'Solid media', 'Stool sample')), 
         treatment = ifelse(treatment == 'cult_70', 'Ethanol shock (70%)', 
                            ifelse(treatment == 'cult_100', 'Ethanol shock (100%)', 
                                   ifelse(treatment == 'cult_T30', 'Heat shock (30 min)', 
                                          ifelse(treatment == 'cult_T60', 'Heat shock (60 min)', 'Stool sample')))))

filter(otu_cultivation, #treatment == c('cult_100', 'cult_70'), 
       value > 0, Phylum != 'Other') %>%
  group_by(treatment, cultivation_media, bile_acids, cultivation_day, Phylum) %>% 
  reframe(sumPA = sum(PA)) %>%  
  ggplot(aes(x = as.factor(cultivation_day), y = sumPA, color = treatment, shape = bile_acids)) +
  geom_jitter(size = 4, width = 0.3) +
  #scale_y_log10() +
  facet_grid(Phylum ~ cultivation_media, scales = 'free', space = 'free_x') +
  labs(x = 'Cultivation day', y = '# OTUs', shape = '', color = '')
ggsave('etoh_t_comparison/plots/number_otus_treatment.png')

cultivationPA <- filter(otu_cultivation, value > 0, Phylum != 'Other') %>%  
  group_by(treatment, cultivation_media, bile_acids, cultivation_day, Phylum) %>% 
  reframe(sumPA = sum(PA))   

left_join(filter(cultivationPA, cultivation_media != 'Stool sample'), 
          filter(cultivationPA, cultivation_media == 'Stool sample'), by = 'Phylum') %>% 
  mutate(ratio_OTUs = sumPA.x/sumPA.y) %>%  
  ggplot(aes(x = as.factor(cultivation_day.x), y = ratio_OTUs, color = treatment.x, shape = bile_acids.x)) +
  geom_jitter(size = 4, width = 0.3) +
  #scale_y_log10() +
  facet_grid(Phylum ~ cultivation_media.x, scales = 'free', space = 'free_x') +
  labs(x = 'Cultivation day', y = 'Ratio OTUs in stool sample to OTUs after treatment and culture', shape = '', color = '')
ggsave('etoh_t_comparison/plots/ratio_otus.png')

filter(otu_cultivation, #treatment == c('cult_100', 'cult_70'), 
       value > 0, Phylum != 'Other') %>%
 # group_by(treatment, cultivation_media, bile_acids, cultivation_day, Phylum) %>% 
#  reframe(sumPA = sum(PA)) %>%  
  ggplot(aes(x = as.factor(cultivation_day), y = rel_abund, fill = treatment)) +
  geom_boxplot() +
  scale_y_log10() +
  #scale_x_continuous(breaks = c(1, 2, 5, 7,14)) +
  facet_grid(Phylum~ cultivation_media + bile_acids , scales = 'free_x', space = 'free_x') +
  labs(x = 'Cultivation day', y = 'Relative abundance', shape = '', fill = '')
ggsave('etoh_t_comparison/plots/relabund_otus_treatment.png')

# Only Bacillota 
filter(otu_cultivation, Phylum == 'Bacillota', value > 0) %>% 
  group_by(treatment, cultivation_media, bile_acids, cultivation_day, Family) %>% 
  reframe(sumPA = sum(PA)) %>%  
  ggplot(aes(x = as.factor(cultivation_day), y = sumPA, color = treatment, shape = bile_acids)) +
  geom_jitter(size = 4, width = 0.3) +
  #scale_y_log10() +
  facet_grid(cultivation_media ~ Family, scales = 'free', space = 'free_y') +
  labs(x = 'Cultivation day', y = '# OTUs', shape = '', color = '')


# Percentages of phyla in
# stool sample
otu_long %>% filter(Group %in% c('UM', 'TZ') & !is.na(Phylum)) %>%
  group_by(Phylum) %>%
  reframe(abund = sum(rel_abund)/2) 

# Heat shock
heat30 <- otu_cultivation %>%
  filter(treatment == 'Heat shock (30 min)') %>%
  group_by(cultivation_day, cultivation_media, bile_acids, Phylum) %>%
  reframe(abund = sum(rel_abund)/n_distinct(Group))

heat60 <- otu_cultivation %>%
  filter(treatment == 'Heat shock (60 min)') %>%
  group_by(cultivation_day, cultivation_media, bile_acids, Phylum) %>%
  reframe(abund = sum(rel_abund)/n_distinct(Group))

# Ethanol shock 
et100 <- otu_cultivation %>%
  filter(treatment == 'Ethanol shock (100%)') %>%
  group_by(cultivation_day, cultivation_media, bile_acids, Phylum) %>%
  reframe(abund = sum(rel_abund)/n_distinct(Group))

et70 <- otu_cultivation %>%
  filter(treatment == 'Ethanol shock (70%)') %>%
  group_by(cultivation_day, cultivation_media, bile_acids, Phylum) %>%
  reframe(abund = sum(rel_abund)/n_distinct(Group))

heat30
heat60
et100
et70

# A plot more for me 
otu_cultivation %>%
  group_by(treatment, cultivation_day, cultivation_media, bile_acids, Phylum) %>%
  reframe(abund = sum(rel_abund)/n_distinct(Group)) %>%
  ggplot(aes(x = cultivation_day, y = abund, color = Phylum, shape = cultivation_media)) +
  geom_point(size = 3) +
  geom_line(mapping = aes(linetype = cultivation_media)) +
  facet_grid(treatment ~ bile_acids * cultivation_media)
ggsave('etoh_t_comparison/plots/abundance_phylum_time.png', dpi = 600)

# 
otu_cultivation %>%
  filter(Phylum == 'Bacillota') %>%
  group_by(treatment, cultivation_day, cultivation_media, bile_acids) %>%
  reframe(abund = sum(rel_abund)/n_distinct(Group)) %>%
  ggplot(aes(x = cultivation_day, y = abund, color = treatment, shape = bile_acids)) +
  geom_point(size = 3, position = position_dodge(width = .3)) +
  geom_line(mapping = aes(linetype = bile_acids), linewidth = 1, position = position_dodge(width = .3)) +
  facet_grid(~ cultivation_media) +
  labs(x = 'Cultivation day', y = 'Relative abundance', color = 'Treatment', shape = 'Bile acids', linetype = 'Bile acids', 
       subtitle = 'Bacillota')
ggsave('etoh_t_comparison/plots/bacillota_cultivation.png', dpi=600)

# 
# #################### HEATMAP with DENDROGRAM ##################
# otu_heat <- otu_long %>% 
#   group_by(treatment, Phylum) %>%
#   summarise(abundmean = mean(rel_abund)) %>%
#   ungroup() %>%
#   filter(!is.na(treatment))
# 
# ggplot(otu_heat, aes(x = Phylum, y = treatment, fill = log10(abundmean))) +
#   geom_tile() +
#   scale_fill_gradient2(low = "yellow", mid = 'red' , high = "darkgreen", midpoint = -3, na.value = 'white') +
#   labs(x = '', y = '') +
#   scale_y_discrete(labels=c('Ethanol shock (100%) + cultivation', 'Ethanol shock (70%) + cultivation','Heat shock (30 min) + cultivation', 'Heat shock (60 min) + cultivation',
#                             'Ethanol shock + EMA 10 min', 'Heat shock + EMA 10 min', 'Ethanol shock + EMA 15 min', 'Heat shock + EMA 15 min', 'Ethanol shock + EMA 5 min',
#                             'Heat shock + EMA 5 min', 'Microbiota', 'Ethanol shock + wash',  'Heat shock + wash')) 
# 
# ggsave('etoh_t_comparison/plots/treatment_phylum.png',dpi= 600)


# Alpha diversity 
richness <- estimateR(otutab)
evenness <- diversity(otutab)/log(specnumber(otutab))

alpha_meta <- as_tibble(as.list(evenness)) %>% pivot_longer(names_to = 'Group', values_to = 'evenness', cols = everything()) %>%
  left_join(t(richness) %>% as.data.frame() %>% rownames_to_column('Group'), by='Group') %>%
  left_join(metadata, by='Group') %>%
  mutate(bile_acids = ifelse(bile_acids == 'yes', 'With bile acids', 'Without bile acids'), 
         cultivation_media = ifelse(cultivation_media == 'liquid', 'Liquid media', 'Solid media'), 
         treatment = ifelse(treatment == 'cult_70', 'Ethanol shock (70%)', 
                            ifelse(treatment == 'cult_100', 'Ethanol shock (100%)', 
                                   ifelse(treatment == 'cult_T30', 'Heat shock (30 min)', 'Heat shock (60 min)'))))

even <- alpha_meta %>%
  filter(cultivation %in% c('cultivationSolid', 'cultivationLiquid')) %>%
  group_by(treatment, cultivation_media, cultivation_day, bile_acids) %>%
  reframe(evenness = mean(evenness, na.rm = TRUE)) %>%
  ggplot(aes(x = cultivation_day, y = evenness, color = treatment)) +
  geom_line(linewidth = 1) +
  #scale_fill_manual(values = col_phylum) +
  facet_grid(cultivation_media ~ bile_acids) +
  labs(x = 'Cultivation day', y = 'Evenness', color = 'Treatment')

obs <- alpha_meta %>%
  filter(cultivation %in% c('cultivationSolid', 'cultivationLiquid')) %>%
  group_by(treatment, cultivation_media, cultivation_day, bile_acids) %>%
  reframe(observed = mean(S.obs, na.rm = TRUE)) %>%
  ggplot(aes(x = cultivation_day, y = observed, color = treatment)) +
  geom_line(linewidth = 1) +
  #scale_fill_manual(values = col_phylum) +
  facet_grid(cultivation_media ~ bile_acids, scales = 'free_y') +
  labs(x = 'Cultivation day', y = '# OTUs', color = 'Treatment')

ggarrange(even + labs(tag = 'A'), obs + labs(tag = 'B'), 
          common.legend = TRUE, legend = 'bottom')
ggsave('etoh_t_comparison/plots/alpha_cultivation.png', dpi=600)


# Beta diversity 
bray <- vegdist(otutab, method = 'bray')

nmds <- metaMDS(bray)
nmds_positions <-
  as.data.frame(scores(nmds, display='sites')) %>%
  rownames_to_column('Group') %>%
  left_join(metadata, by = 'Group') 

stool <- filter(nmds_positions, Group %in% c('UM', 'TZ'))

ord_cultivation <- nmds_positions %>%
  filter(cultivation %in% c('cultivationSolid', 'cultivationLiquid')) %>%
  rbind(stool)


ord_person <- ord_cultivation %>%
  ggplot(aes(x=NMDS1, y=NMDS2, color=person)) +
  geom_point(size=3) +
  #scale_color_manual(values = c('darkgreen', 'darkblue')) +
  labs(x = '', y='', color = 'Individual', tag= 'A')

ord_media <- ord_cultivation %>%
  ggplot(aes(x=NMDS1, y=NMDS2, color=cultivation_media)) +
  geom_point(size=3) +
  #scale_color_manual(values = c('darkred', 'skyblue', 'gold3')) +
  labs(x = '', y='', color = 'Cultivation media', tag = 'C')

ord_bile <-  ord_cultivation %>%
  ggplot(aes(x=NMDS1, y=NMDS2, color=bile_acids)) +
  geom_point(size=3) +
  #scale_color_manual(values = c('orchid4', 'green4')) +
  labs(x = '', y='', color = 'Bile acids', #tag = 'D'
       )

ord_tretament <-  ord_cultivation %>%
  ggplot(aes(x=NMDS1, y=NMDS2, color=treatment)) +
  geom_point(size=3) +
  #scale_color_manual(values = c('green', 'red', 'dodgerblue', 'yellow2', 'hotpink3')) +
  labs(x = '', y='', color = 'Treatment', tag = 'B')

ord_day <-  ord_cultivation %>%
  ggplot(aes(x=NMDS1, y=NMDS2, color=as.factor(cultivation_day))) +
  geom_point(size=3) +
  labs(x = '', y='', color = 'Cultivation day', tag = 'E')

ggarrange(ord_person, ord_tretament, ord_media, ord_bile, ord_day, common.legend = FALSE)
ggsave('etoh_t_comparison/plots/beta_cultivation.png', dpi=600)

# # PCA plot
# otu_long  <- pivot_longer(as.data.frame(otutab) %>%  
#                             rownames_to_column('sample'), -sample)
# 
# pca <- vegan::rda(otutab)
# 
# biplot(pca,
#        display = c("sites", 
#                    "species"),
#        type = c("text",
#                 "points"))

## Is there any correlation between the CFU on plates and the observed number of OTUs? 
unique(alpha_meta$description)
unique(cfu$description)

alpha_cfu <- alpha_meta %>% left_join(cfu, by = c('description' ,'experiment', 'shock', 
                                                  'treatment', 'cultivation_day', 'person', 
                                                  'cultivation', 'cultivation_media', 'bile_acids'))

et_cor <- cor.test(filter(alpha_cfu, shock == 'Ethanol_shock')$S.obs, filter(alpha_cfu, shock == 'Ethanol_shock')$CFU, use = 'complete.obs', method = 'pearson')
h_cor <-  cor.test(filter(alpha_cfu, shock == 'Heat_shock')$S.obs, filter(alpha_cfu, shock == 'Heat_shock')$CFU, use = 'complete.obs', method = 'pearson')

alpha_cfu %>%
  filter(!is.na(shock)) %>%  
  mutate(shock = ifelse(shock == 'Ethanol_shock', 'Ethanol shock', 
                        ifelse(shock == 'Heat_shock', 'Heat shock', 'Stool sample'))) %>%
  ggplot(aes(x = CFU, y = S.obs, color = shock)) +
  geom_jitter(size = 3, width=5) +
  scale_color_manual(values = col_shock) +
  annotate('text', y = 200, x = 200, label = paste('Ethanol shock', '\n','Correlation:', round(et_cor$estimate, digits = 2), '\n', 
                                                   'p =', scientific(et_cor$p.value, digits =1))) +
  annotate('text', y = 250, x = 200, label = paste('Heat shock', '\n','Correlation:', round(h_cor$estimate, digits = 2), '\n', 
                                                   'p =', scientific(h_cor$p.value, digits =1))) +
  labs(x = 'CFU/ml', y='#OTUs', color = '') +
  geom_smooth(method = 'lm')
ggsave('etoh_t_comparison/plots/corr_observed_CFU_smooth.png', dpi=600) 



############ EMA treatment ############# 
# For EMA treatments sepaeratly 
conc_ema = metadata %>% 
  filter(cultivation == 'no') %>%
  filter(!(Group %in% c('UM', 'TZ'))) %>%
  group_by(treatment) %>%
  mutate(meanConc = mean(DNAconcentration)) %>%
  ungroup()

conc_ema$treatment <- factor(conc_ema$treatment, levels = c('none', 'wash_T', 'EMA5_T', 'EMA10_T', 'EMA15_T', 
                                                            'wash_Et', 'EMA5_Et', 'EMA10_Et', 'EMA15_Et')) 

conc_ema %>% ggplot(aes(x=DNAconcentration, y=treatment, fill = treatment, color = treatment)) +
  geom_point(size=3) +
  geom_boxplot(alpha = .3, outlier.shape = 8, outlier.color = 'black', outlier.size = 5) +
  scale_y_discrete(labels=c('Stool sample', 'Heat shock + wash', 'Heat shock + EMA 5 min', 'Heat shock + EMA 10 min', 'Heat shock + EMA 15 min', 
                            'Ethanol shock + wash', 'Ethanol shock + EMA 5 min','Ethanol shock + EMA 10 min', 'Ethanol shock + EMA 15 min')) +
  labs(x= 'DNA concentration (ng/μl)', y='') +
  theme(legend.position="none")
ggsave('etoh_t_comparison/plots/DNAconcentration_ema.png', dpi=600)


metadata %>% 
  filter(cultivation == 'no') %>%
  filter(!(Group %in% c('UM', 'TZ'))) %>%  
  group_by(treatment) %>%  
  reframe(n = n_distinct(Group))
#
otu_ema <- filter(otu_long, !is.na(treatmentEMA)) %>%
  select(Group, description, shock, treatmentEMA, aliquote, name, value, rel_abund, PA, Phylum, Class, Order, Family, Genus)

unique(otu_ema$Group)

# relative abundance plot 
otu_ema_rel <- filter(otu_ema, shock == 'Non-treated') %>%  
  group_by(name, Phylum, Class, Order, Family, Genus) %>%  
  reframe(rel_abund = mean(rel_abund)) %>%
  full_join(filter(otu_ema, shock != 'Non-treated'), by = c('name', 'Phylum', 'Class', 'Order', 'Family', 'Genus'), relationship = "many-to-many") %>%
  mutate(shock = case_when(shock == 'Ethanol_shock' ~ 'Ethanol shock', 
                           shock == 'Heat_shock' ~ 'Heat shock', 
                           shock == 'Only-water' ~ 'Only water'))

otu_ema_rel %>%
  filter(rel_abund.x != 0 & rel_abund.y != 0) %>%
  mutate(ratio = rel_abund.x/rel_abund.y) %>%
  group_by(Phylum, shock, treatmentEMA) %>%
  reframe(mean_ratio = mean(ratio)) %>%
  ggplot(aes(x = mean_ratio, y = Phylum, color = as.factor(treatmentEMA))) +
  #geom_boxplot() +
  scale_x_log10() +
  geom_jitter(size = 4, width = .2) +
  geom_vline(xintercept = 1) +
  facet_grid(~ shock) +
  labs(x = 'Ratio of relative abundance in stool sample & treated sample', y = '', color = 'Illumination time \n for EMA treatment')
ggsave('etoh_t_comparison/plots/efficiency_EMA_treatment_class.png', dpi=600)


otu_ema2 <- otu_ema %>% 
  mutate(treatmentEMA = case_when(treatmentEMA == '5' ~'5 min', 
                                  treatmentEMA == '10' ~'10 min', 
                                  treatmentEMA == '15' ~'15 min'), 
         treatmentEMA = ifelse(shock == 'Only-water', 'Only water', 
                               ifelse(shock == 'Non-treated', 'Non-treated', treatmentEMA))) 
otu_ema2$treatmentEMA <- factor(otu_ema2$treatmentEMA, 
                                levels = c('Non-treated', 'Only water', '5 min', 
                                            '10 min', '15 min'))

otu_ema2 %>% 
  filter(Phylum != 'Other') %>% 
  ggplot(aes(x = treatmentEMA, y = rel_abund, fill = shock)) +
  geom_boxplot() +
  scale_y_log10() +
  facet_wrap(~ Phylum, scales = 'free_y', ncol = 3) +
  labs(y = 'Relative abundance [log10]', y = '', fill = '')
ggsave('etoh_t_comparison/plots/rel_abundEMA.png')

# Beta diversity (to check how is the reproducibility of treatment)
ema_bray <- filter(nmds_positions, !is.na(treatmentEMA))

ema_bray %>%
  ggplot(aes(x=NMDS1, y=NMDS2, color=as.factor(treatmentEMA), shape = shock)) +
  geom_point(size=4) +
  #scale_color_manual(values = c('darkred', 'skyblue', 'gold3')) +
  labs(x = '', y='', color = 'Time of light \n incubation with EMA', shape = '')
ggsave('etoh_t_comparison/plots/bray_ema.png', dpi=600)




