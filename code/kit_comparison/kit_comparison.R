# Analysis of V3V4 16S rRNA data for comparison of kits for DNA isolations of human stool samples 
# OTUs were constructed in mothur, code availabile HOPC/bin/mothur.script 

library(dplyr)
library(tidyr)
library(tibble)
library(vegan)
library(stringr)
library(readr)
library(ggplot2)
library(scales)
library(ggpubr)

set.seed(96)
theme_set(theme_bw())

# Exploration
shared = read.table('kit_comparison/data/mothur/stability.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.opti_mcc.shared', sep = '\t', header = TRUE) %>%
  select(Group, starts_with('Otu')) %>%
  pivot_longer(-Group)

# Distribution of reads per sample 
shared %>%
  group_by(Group) %>%
  summarize(ReadsPerSample=sum(value)) %>%
  ggplot(aes(x=ReadsPerSample)) +
  geom_histogram() +
  labs(x = 'Reads per sample', y = 'Number of samples')
ggsave('kit_comparison/plots/reads_per_sample.png', dpi=600)

# How min/max/mean/sum reads per sample/all samples
info1 = shared %>%
  group_by(Group) %>%
  summarize(ReadsPerSample=sum(value)) %>% 
  # remove negative controls
  filter(!(Group %in% c('NCD', 'NCF', 'NCN', 'NCP')))
min(info1$ReadsPerSample) # 12 542
max(info1$ReadsPerSample) # 153 530
mean(info1$ReadsPerSample) # 90515.81
median(info1$ReadsPerSample) # 90499.5
sum(info1$ReadsPerSample) # 3258569

# From the plot above I see that 50 000 will be a good rarefaction depth, as only 5 samples are below this 
# ( 4 negative controls + 1 sample)

# Distribution of OTUs
shared %>%
  group_by(name) %>%
  summarize(OTUabundance=sum(value)) %>%
  ggplot(aes(x=OTUabundance)) +
  geom_histogram(breaks=seq(0, 250, by =1)) +
  labs(x = 'Abundance of OTU', y= 'Number of OTUs')
ggsave('kit_comparison/plots/reads_per_OTU.png', dpi=600)

# Number of OTUs
length(unique(shared$name)) # I have 14 045 OTUs

# How min/max/mean reads per OTU
info2 = shared %>%
  group_by(name) %>%
  summarize(OTUabundance=sum(value))
min(info2$OTUabundance) # 1
max(info2$OTUabundance) # 208785
mean(info2$OTUabundance) # 232.0283
median(info2$OTUabundance) # 1
sum(info2$OTUabundance) # 3258837 = Total number of redas in the dataset before removing anything! 

reads_per_OTU = shared %>%
  group_by(name) %>%
  summarize(OTUabundance=sum(value))

# How many OTUs have less than 5 reads sum; sumes up TRUE (as this is 1 in computer code)
sum(reads_per_OTU$OTUabundance < 5) # 12776
# How many singletons ? 
reads_per_OTU %>%
  filter(OTUabundance == 1) %>%
  reframe(n_otus = n_distinct(name),  # 11 312 80% OTUs are singletons, but their abundance 0.35%
          n_reads = (n_otus/3258837)*100) 

# How many reads do they contain
reads_per_OTU %>%
  filter(OTUabundance < 5) %>%
  summarise(sum(OTUabundance)) # 15061 = 0.46% of all reads

sum(reads_per_OTU$OTUabundance < 10) # 13120
reads_per_OTU %>%
  filter(OTUabundance < 10) %>%
  summarise(sum(OTUabundance)) # 17257 = 0.53% of all reads

# Min number of reads in sample and min sequences per OTU as determined in the exploration analysis
reads_per_sample = 50000
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
saveRDS(otutab, 'kit_comparison/data/r_data/otutab.RDS')

# Extract OTUs that are present rarefied table 
otu_names = as.data.frame(otutab) %>% colnames() 

# Import taxonomy table
taxtab = read.table('kit_comparison/data/mothur/stability.trim.contigs.good.unique.good.filter.unique.precluster.denovo.vsearch.pick.opti_mcc.0.03.cons.taxonomy', sep = '\t', header = TRUE) %>%
  filter(OTU %in% otu_names) %>%
  select(name = "OTU", taxonomy = "Taxonomy") %>%
  mutate(taxonomy = str_replace_all(taxonomy, "\\\\|\\\"|\\(\\d+\\)", ""),
         taxonomy = str_replace(taxonomy, ";$", "")) %>%
  separate(taxonomy, into=c("Domain", "Phylum", "Class", "Order", "Family", "Genus"),
           sep=";") 
saveRDS(taxtab, 'kit_comparison/data/r_data/taxtab.RDS')

# Import metadata
metadata = as_tibble(read.csv('kit_comparison/data/metadata.csv', sep=',')) %>%
  filter(Group %in% rownames(otutab))

# remove unnecessary 
rm(info1)
rm(info2)
rm(min_seqs_per_otu)
rm(otu_names)
rm(reads_per_sample)
rm(reads_per_OTU)
rm(shared)

# DNA concentrations 
col_kit = c("#D22B2B", "#2E9FDF", "#E7B800", '#27ae60')

conc = as_tibble(read.csv('kit_comparison/data/concentration.csv', sep=',')) %>%
  filter(Organism != 'Negative control') %>%
  mutate(Concentration = if_else(Concentration < 0, 0, Concentration))# write 0 if the concentration is under the detection limit

anova_res <- aov(Concentration ~ Protocol, data = conc) 
summary(anova_res)

res <- pairwise.t.test(conc$Concentration, conc$Protocol, paired = TRUE, p.adjust.method = 'BH') 

res_t <- as.data.frame(res$p.value) %>%
  rownames_to_column('Protocol1') %>%
  pivot_longer(values_to = 'pvalue', names_to = 'Protocol2', cols = c(2:4)) %>%
  mutate(Protocol = paste(Protocol1,'&',Protocol2)) %>%
  filter(Protocol1 != Protocol2) %>%
  na.omit()

conc_plot <- ggplot(conc, mapping = aes(x=Protocol, y=Concentration, fill=Protocol)) + 
  geom_boxplot() +
  stat_compare_means(method = 't.test', paired = TRUE, mapping = aes(label = paste0("p = ", after_stat(p.format))), 
                     comparisons = list(c('FastDNA', 'PowerFecal'), c('PowerFecal', 'NucleoSpin'), c('PowerFecal', 'DNAeasy')
                                        # ,c('FastDNA', 'NucleoSpin'), c('FastDNA', 'DNAeasy'), c('NucleoSpin', 'DNAeasy')
                                        )) +
  scale_fill_manual(values = col_kit) +
  labs(x='', y="DNA concentration (ng/µl)", color="") +
  theme(legend.position = 'none') 
conc_plot
ggsave('kit_comparison/plots/DNA_concentrations.png', dpi=600)

# Alpha diversity 
# Calculated with relative abundances! 
richness = estimateR(otutab) # observed richness and Chao1
evenness = diversity(otutab)/log(specnumber(otutab)) # evenness index
shannon = diversity(otutab, index = 'shannon')

# Join all calculations and metadata
alpha_meta = as_tibble(as.list(evenness)) %>% pivot_longer(names_to = 'Group', values_to = 'evenness', cols = starts_with('H')) %>%
  left_join(t(richness) %>% as.data.frame() %>% rownames_to_column('Group'), by='Group') %>%
  left_join(as_tibble(as.list(shannon)) %>% pivot_longer(names_to = 'Group', values_to = 'shannon', cols = starts_with('H'))) %>%
  left_join(metadata, by='Group') %>%
  na.omit()

evenness_plot <- ggplot(alpha_meta, mapping = aes(x=kit, y=evenness, color=kit)) + 
  geom_point(size = 3) + 
  stat_compare_means(mapping = aes(label = after_stat(p.signif)), method = "t.test", paired = FALSE,
                     comparisons=list(c('FastDNA', 'Dneasy'), c('PowerFecal', 'Dneasy'), c('NucleoSpin', 'Dneasy')
                                      #c("NucleoSpin", "FastDNA"), c("NucleoSpin", "PowerFecal"), c("FastDNA", "PowerFecal"), 
                                      )) +
  scale_color_manual(values = col_kit) +
  labs(x='', y="Evenness", color="") +
  theme(legend.position = 'none') 
evenness_plot

# Chao1 
chao_plot <- ggplot(alpha_meta, mapping = aes(x=kit, y=S.chao1, color=kit)) + 
  geom_point() + 
  # Nothing is significant
  # stat_compare_means(mapping = aes(label = after_stat(p.signif)), method = "t.test", paired=F,
  #                    comparisons=list(c("NucleoSpin", "FastDNA"), c("NucleoSpin", "PowerFecal"), c("FastDNA", "PowerFecal"), 
  #                                     c('FastDNA', 'Dneasy'), c('PowerFecal', 'Dneasy'), c('NucleoSpin', 'Dneasy'))) +
  scale_color_manual(values = col_kit) +
  labs(x='', y="Chao1", color="") +
  theme(legend.position = 'none') 
chao_plot

# Observed 
observed_plot <- ggplot(alpha_meta, mapping = aes(x=kit, y=S.obs, color=kit)) + 
  geom_point() + 
  # Nothing is significant
  # stat_compare_means(mapping = aes(label = after_stat(p.signif)), method = "t.test", paired=F,
  #                    comparisons=list(c("NucleoSpin", "FastDNA"), c("NucleoSpin", "PowerFecal"), c("FastDNA", "PowerFecal"), 
  #                                     c('FastDNA', 'Dneasy'), c('PowerFecal', 'Dneasy'), c('NucleoSpin', 'Dneasy'))) +
  scale_color_manual(values = col_kit) +
  labs(x='', y="Observed OTUs", color="") +
  theme(legend.position = 'none') 
observed_plot

# Shannon
alpha_shannon <- pairwise.t.test(alpha_meta$shannon, alpha_meta$kit, paired = FALSE, p.adjust.method = 'BH')
p_shannon <- pairwise.t.test(conc$Concentration, conc$Protocol, paired = TRUE, p.adjust.method = 'BH') 

shannon <- ggplot(alpha_meta, mapping = aes(x=kit, y=shannon, color=individual)) + 
  geom_point(size = 3) + 
  stat_compare_means(mapping = aes(label = after_stat(p.signif)), method = "t.test", paired=F,
                     comparisons=list(#c("NucleoSpin", "FastDNA"), c("NucleoSpin", "PowerFecal"), c("FastDNA", "PowerFecal"),
                                      c('FastDNA', 'Dneasy'), c('PowerFecal', 'Dneasy'), c('NucleoSpin', 'Dneasy'))) +
  #scale_color_manual(values = col_kit) +
  labs(x='', y="Shannon", color="") +
  theme(legend.position = 'bottom') 
shannon

# Beta diversity 
dist <- vegdist(otutab, method = 'bray')

nmds <- metaMDS(dist)
nmds_positions <-
  as.data.frame(scores(nmds, display='sites')) %>%
  rownames_to_column('Group') %>%
  left_join(metadata, by = 'Group') %>%
  na.omit()

permanova <- adonis2(dist ~ kit, data = metadata, method = 'bray', permutations = 999)

kits <- ggplot(nmds_positions, aes(x=NMDS1, y=NMDS2, color=kit, shape = individual)) +
  geom_point(size=4) +
  annotate('text', x = 0.5, y = 0.6, label = paste('PERMANOVA, p = 0.002')) +
  scale_color_manual(values = col_kit) +
  labs(x='', y='', color='Protocol', shape = 'Individual')
kits
ggsave('kit_comparison/plots/nmds_bray.png', dpi=600)


# Taxonomy 
rel_abund <- otutab %>% 
  as.data.frame() %>% 
  rownames_to_column('Group') %>% 
  pivot_longer(starts_with('Otu')) %>%
  left_join(metadata, by = 'Group') %>% 
  left_join(taxtab, by = 'name') %>% 
  group_by(Group) %>% 
  mutate(rel_abund = value/sum(value)) %>% 
  ungroup() %>%  
  group_by(kit, Phylum) %>% 
  summarise(sum = sum(rel_abund)/n_distinct(Group) * 100) %>% 
  ungroup() %>% 
  filter(sum > 0.01) %>% 
  ggplot(aes(x = kit, y = sum, fill = Phylum)) +
  geom_col() +
  labs(x = '', y = 'Relative abundance [%]')
rel_abund

ggarrange(conc_plot + labs(tag = 'A'), shannon + labs(tag = 'B'),  
          kits + labs(tag = 'C'), rel_abund + labs(tag = 'D'),
          nrow = 2, ncol = 2)
ggsave('kit_comparison/plots/kit_comparison_v2.png', dpi = 600) 



# Percentage of shared OTUs based on rarefaction /sequencing depth 
sequencing_depths <- c(10000, 20000, 30000, 40000, 50000, 60000, 70000, 80000)
meta_pre <- as_tibble(read.csv('kit_comparison/data/metadata.csv', sep=',')) 

n_otus_shared = data.frame()
for(i in sequencing_depths) {
  shared_pre <- shared %>%
    group_by(Group) %>%
    mutate(sum_sample = sum(value)) %>%
    ungroup() %>%
    filter(sum_sample > i) %>%
    filter(!(Group %in% c('NCD', 'NCF', 'NCN', 'NCP'))) %>%
    select(-sum_sample)
  
  meta <- meta_pre %>%
    filter(Group %in% shared_pre$Group)

  otutab_pre <- shared_pre %>%
    pivot_wider(names_from = 'name', values_from = 'value', values_fill = 0) %>%
    column_to_rownames('Group')
  
  rarefied <- vegan::rrarefy(otutab_pre, sample = i)
  
  long <- as.data.frame(rarefied) %>%
    rownames_to_column('Group') %>%
    pivot_longer(names_to = 'name', values_to = 'value', cols = starts_with('Otu')) %>%
    left_join(meta, by = 'Group')
  
  n_all <- long %>%
    filter(value > 0) %>%
    group_by(kit, individual) %>%
    summarise(n_all = n_distinct(name), .groups = 'drop')
  
  n_shared <- long %>% 
    mutate(value = ifelse(value > 0, 1, 0)) %>%
    filter(value > 0) %>% 
    group_by(name, kit, individual) %>%
    summarise(value = ifelse(sum(value) == 3, 1, 0), .groups = 'drop') %>%
    group_by(kit, individual) %>%
    summarise(n_otus = sum(value), .groups = 'drop') %>%
    left_join(n_all, by = c('kit', 'individual')) %>%
    mutate(seq_depth = i, percent = (n_otus/n_all) *100) %>%
    filter(n_otus > 0)
  
  n_otus_shared <- rbind(n_otus_shared, n_shared)
  
}

rarefaction <- ggplot(n_otus_shared, aes(x = seq_depth, y = n_otus, color = kit, shape = individual)) +
  geom_jitter(size = 3, width = 1000) +
  scale_color_manual(values = col_kit) +
  labs(x = 'Rarefaction depth', y = '#OTUs shared', color = 'Protocol', shape = '')


# What was in the negative control samples ?
metadata = as_tibble(read.csv('kit_comparison/data/metadata.csv', sep=','))

nc <- shared %>%  
  left_join(metadata, by = 'Group') %>% 
  left_join(taxtab, by = 'name') %>% 
  filter(group == 'negative_control', value > 0)


