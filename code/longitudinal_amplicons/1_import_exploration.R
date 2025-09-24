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
theme_set(theme_bw())

shared <- read.table('../longitudinal_amplicons/data/mothur/final.opti_mcc.shared', sep = '\t', header = TRUE) %>%
  select(Group, starts_with('Otu')) %>%
  pivot_longer(-Group)

# Distribution of reads per sample 
shared %>%
  group_by(Group) %>%
  summarize(ReadsPerSample=sum(value)) %>%
  ggplot(aes(x=ReadsPerSample)) +
  geom_histogram() +
  scale_x_continuous(breaks = seq(0, 800000, by=100000)) +
  scale_y_continuous(breaks = seq(0,40, by=5)) +
  labs(x = 'Reads per sample', y = 'Number of samples')
ggsave('thesis/plots/reads_per_sample.png', dpi=600)

# How min/max/mean/sum reads per sample/all samples
info1 = shared %>%
  group_by(Group) %>%
  summarize(ReadsPerSample=sum(value)) %>% 
  filter(Group != 'SNC') %>% filter(Group != 'MNC')
min(info1$ReadsPerSample) #61388
max(info1$ReadsPerSample) #884743
mean(info1$ReadsPerSample) #270662.7
median(info1$ReadsPerSample) # 237272.5
sum(info1$ReadsPerSample) #62793747

# Distribution of OTUs
shared %>%
  group_by(name) %>%
  summarize(OTUabundance=sum(value)) %>%
  ggplot(aes(x=OTUabundance)) +
  geom_histogram(breaks=seq(0, 20, by =1)) +
  labs(x = 'Number of reads per OTU', y= 'Number of OTUs')
ggsave('thesis/plots/reads_per_OTU.png', dpi=600)

# How min/max/mean reads per OTU
info2 = shared %>%
  group_by(name) %>%
  summarize(OTUabundance=sum(value))
min(info2$OTUabundance) # 1
max(info2$OTUabundance) # 14087376
mean(info2$OTUabundance) # 310.0056
median(info2$OTUabundance) # 1
sum(info2$OTUabundance) #62797212

# Total number of redas in the dataset before removing anything! 
sum(shared$value) #62226327

reads_per_OTU = shared %>%
  group_by(name) %>%
  summarize(OTUabundance=sum(value))

# How many OTUs have less thaN 10 reads
length(reads_per_OTU$OTUabundance < 10) # 202568
# How many reads do they contain
sum(reads_per_OTU$OTUabundance < 10) # 191964
# what percentage is this from all reads ?
sum(reads_per_OTU$OTUabundance < 10)/sum(reads_per_OTU$OTUabundance)*100 # 0.3

# We loose 202568 OTUs, but this represents only 0.3 % of all reads 

sum(reads_per_OTU$OTUabundance < 62)/sum(reads_per_OTU$OTUabundance) *100 # 31% readov. 
# Removing rare OTUs, what is happening with the data 

# How many reads / OTUs I loose removing xyxy rare OTUs?
# 0%, 100%, 10%, 1%, 0.1%, 0.01%, 0.001%, 0.0001%, 0.00001%, 0,000001%
percent_removed = c(0, 1, 0.1, 0.01, 0.001, 0.0001, 0.00001, 0.000001, 0.0000001, 0.00000001)
results = data.frame()

for (i in percent_removed) {
  number_reads = sum(shared$value) * i
  
  shared1 = shared %>%
    group_by(name) %>%
    summarise(sum_reads = sum(value), .groups = 'drop') %>%
    filter(sum_reads > number_reads)
  
  res_sum = sum(shared1$sum_reads)
  res_distinct = n_distinct(shared1$name)
  
  results = rbind(results, data.frame(percent_removed = i, 
                                      number_reads = res_sum, 
                                      number_otus = res_distinct))
  
}

write_csv(results, 'thesis/data/removed.csv')

ggarrange(ggplot(results, aes(x = percent_removed, y = number_reads)) +
            geom_point(size = 2) +
            scale_x_continuous(labels = function(x) paste0(x, "%")) +
            labs(x = 'Percent of reads removed', y = 'Number of reads remaining'), 
          ggplot(results %>%
                   filter(percent_removed > 1e-7), aes(x = percent_removed , y= number_otus)) +
            geom_point(size = 2) +
            scale_x_continuous(labels = function(x) paste0(x, "%")) +
            labs(x = 'Percent of reads removed', y = 'Number of OTUs remaining'))
ggsave('thesis/plots/otu_reads_removed.png', dpi = 600)

# negative control samples 
shared %>%  filter(Group %in% c('MNC', 'SNC')) %>% 
  group_by(Group) %>% 
  reframe(n_reads = sum(value))

n <- shared %>%  filter(Group %in% c('MNC', 'SNC')) %>% 
  group_by(Group, name) %>% 
  reframe(n_otu = sum(value)) %>%  
  filter(n_otu > 0)


###
## Rarefaction 
###
# 1. Difference between exact and empirical rarefaction 

# If we choose a sampling depth how many OTUs will be observed at that point?
# A function to subsample data once 
subsample = function(data, sample_size){
  data %>%
    group_by(Group) %>%
    uncount(value) %>%
    sample_n(size=sample_size) %>%
    summarise(noOTU = n_distinct(name))
}

shared2 = shared %>% 
  group_by(Group) %>%
  # Calculate the total number of seqs in each sample and remove samples with less than 150 000
  mutate(total_sample=sum(value)) %>%
  filter(total_sample > 150000) %>%
  group_by(name) %>%
  # Calculate the total number of seqs for each otu and remove otus with less than 10 reads
  mutate(total_otu=sum(value)) %>%
  filter(total_otu > 10) %>%
  ungroup() %>%
  select(-total_sample, -total_otu)

# Do multiple subsamplings
subsamplings = map_dfr(1:100, ~subsample(shared2, 150000), .id = 'iters')

# Calculate from our function the empirical number of OTUs we would get at a choosen sequencing depth
empirical = subsamplings %>% 
  group_by(Group) %>%
  summarise(noOTUs = mean(noOTU))

# Use VEGANS exact number of rarefaction! 
exact = shared2 %>%
  group_by(Group) %>%
  summarise(noOTUs = rarefy(value, 150000))

bind_rows(empirical=empirical, exact=exact, .id='approach') %>%
  ggplot(aes(x=approach, y=noOTUs, group=Group)) +
  geom_line(linewidth = 1) +
  labs(x = '', y = 'Number of OTUs')
ggsave('thesis/plots/empirical_exact_line.png', dpi=600)

inner_join(empirical, exact, by='Group') %>%
  ggplot(aes(x=noOTUs.x, y=noOTUs.y)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, color = 'grey') +
  labs(x='Empirical number of OTUs', y= 'Extact number of OTUs')
ggsave('thesis/plots/empirical_exact_slope.png', dpi=600)

# What is the percent error between empirical and exact
inner_join(empirical, exact, by='Group') %>%
  mutate(error=100*(noOTUs.y-noOTUs.x)/noOTUs.y) %>%
  ggplot(aes(x=error)) +
  geom_density() +
  labs(title = 'What is the percent error between empirical and exact?')
ggsave('thesis/plots/exact_empirical_density.png', dpi=600)

# calculcate the exact mean and sd of error 
inner_join(empirical, exact, by='Group') %>%
  mutate(error=100*(noOTUs.y-noOTUs.x)/noOTUs.y) %>%
  summarise(mean= mean(error), sd=sd(error))
# mean 0.00732 sd 0.143

# 2. Best sampling depth for rarefaction 
# Calculate the sampling coverage for each sample 
sampling_covergae = shared %>%
  group_by(Group) %>%
  summarise(n_seqs=sum(value))

# Plot sampling coverage with different plots to get an idea what my data looks like
sampling_covergae %>% 
  ggplot(aes(x=n_seqs)) +
  geom_histogram(binwidth = 10000) +
  coord_cartesian(xlim=c(0,200000))
ggsave('thesis/plots/sampling_coverage_histo_partial.png', dpi=600)

sampling_covergae %>% 
  ggplot(aes(x=1, y=n_seqs)) +
  geom_jitter() +
  scale_y_log10()
ggsave('thesis/plots/sampling_coverage_jitter.png', dpi=600)

sampling_covergae %>%
  arrange(n_seqs) %>%
  ggplot(aes(x=1:nrow(.), y=n_seqs))+
  geom_line() + 
  coord_cartesian(xlim=c(0,50), ylim=c(0, 200000))
ggsave('thesis/plots/sampling_coverage_line_partial.png', dpi=600)

## Good's coverage 
# Good's coverage = fraction of sequences that appear in an OTU that has been seen more than 1.
coverage_stats = shared %>%
  group_by(Group) %>%
  summarise(n_seqs=sum(value), 
            n_singletons = sum(value==1), 
            goods= 100*(1- n_singletons/n_seqs))

coverage_stats%>%
  ggplot(aes(x=n_seqs, y=goods)) +
  geom_point()+
  labs (x = 'Number of OTUs', y='Goods coverage')
ggsave('thesis/plots/goods_coverage.png', dpi=600)


# 3. Rarefaction effect on diveristy 
# Comparison of the effect of no rarefaction has on ALPHA DIVERSITY METRICS (observed, Chao, Shannon, Simpson) 
# Effect on alpha diveristy measures 
shared %>%
  # Make a randomized otutab, without subsampling !
  uncount(value) %>%                                      
  mutate(name = sample(name)) %>%
  count(Group, name, name='value') %>%
  group_by(Group) %>%
  # Calculate alpha diversity metrics ! 
  summarise(Observed = specnumber(value),               
            Shannon = diversity(value, index = "shannon"), 
            Simpson = diversity(value, index="simpson"), 
            n=sum(value)) %>%
  pivot_longer(cols = c(Observed, Shannon, Simpson), 
               names_to = 'metric') %>%
  ggplot(aes(x=n, y=value)) +
  geom_point() +
  geom_smooth()+
  facet_wrap(~metric, nrow=3, scales='free_y') +
  labs(x = 'Number of OTUs per sample', y= 'Diversity value')
ggsave('thesis/plots/alpha_samplingEffect.png', dpi=600)

# Comparison of the effect of no rarecation / rarefaction / relative abundance /normalization on Bray-Curtis distances
# This matrix has samples that are not statisticly different from eachother
rand = shared %>%
  uncount(value) %>%
  mutate(rand_name = sample(name)) %>%
  select(-name) %>%
  count(Group, rand_name)
# Turn into matrix so it can be used in distance calculations
rand_df = rand %>%
  pivot_wider(names_from="rand_name", values_from="n", values_fill = 0) %>%
  column_to_rownames('Group') %>%
  as.data.frame()

rand_matrix <- rand_df %>%
  as.matrix()

# From rand matrix! 
# Calculate one itteration of distance matrix and avgdist - default is 100 
norare_dist_matrix <- vegdist(rand_matrix, method="bray")
rare_dist_matrix <- avgdist(rand_matrix, dmethod="bray", sample=150000, iterations=9)

norare_dist_tibble <- norare_dist_matrix %>%
  as.matrix() %>%
  as_tibble(rownames="sample") %>%
  pivot_longer(-sample) %>%
  filter(name < sample)

rare_dist_tibble <- rare_dist_matrix %>%
  as.matrix() %>%
  as_tibble(rownames="sample") %>%
  pivot_longer(-sample) %>%
  filter(name < sample)

# Calculate avgdist from relative abundance data
relabund_matrix = rand %>%
  group_by(Group) %>%
  mutate(rel_abund = n/sum(n)) %>%
  ungroup() %>%
  select(-n) %>%
  pivot_wider(names_from = 'rand_name', values_from = 'rel_abund', values_fill = 0) %>%
  as.data.frame() %>%
  column_to_rownames('Group') %>%
  as.matrix()

relabund_dist_matrix = vegdist(relabund_matrix, method='bray')

relabund_dist_tibble <- relabund_dist_matrix %>%
  as.matrix() %>%
  as_tibble(rownames="sample") %>%
  pivot_longer(-sample) %>%
  filter(name < sample)

# Normalization 
# What is the smallest sample (the minimum number of sequences in a sample?)
rand_group_count = rand %>% 
  group_by(Group) %>%
  summarise(n=sum(n))
min_group = min(rand_group_count$n)

# Make a normalizied matrix
# package SRS has built in functions for normalization of ecology data: 
library(SRS)  # https://pubmed.ncbi.nlm.nih.gov/32832266/

# SRS needs a otutab that has OTUs as rows and samples by columns
# Normalized to the min_group number! 
normalized = rand_df %>%
  t() %>%
  as.data.frame() %>%
  SRS(Cmin=min_group)
# Check if all samples have min_group! 
normalized %>% as_tibble(rownames = 'otu') %>%
  pivot_longer(-otu) %>%
  group_by(name) %>%
  summarise(n=sum(value))

normalized_dist_matrix = vegdist(t(normalized), method='bray') 

normalized_dist_tibble <- normalized_dist_matrix %>%
  as.matrix() %>%
  as_tibble(rownames="sample") %>%
  pivot_longer(-sample) %>%
  filter(name < sample) 

comparison <- inner_join(norare_dist_tibble, rare_dist_tibble, by=c("sample", "name")) %>%
  inner_join(., relabund_dist_tibble, by=c('sample', 'name')) %>%
  inner_join(., normalized_dist_tibble, by=c('sample', 'name')) %>%
  select(sample, name, norarefied=value.x, rarefied=value.y, relabund=value.x.x, normalized=value.y.y) %>%
  inner_join(., rand_group_count, by=c("sample" = "Group")) %>%
  inner_join(., rand_group_count, by=c("name" = "Group")) %>%
  mutate(n_diff = abs(n.x-n.y)) %>%
  select(-n.x, -n.y)

comparison %>%
  pivot_longer(cols=c("norarefied", "rarefied", 'relabund', 'normalized'), names_to="type", values_to="dist") %>%
  ggplot(aes(x=n_diff,  y=dist)) +
  geom_point(size=0.25, alpha=0.25) +
  facet_wrap(~type, nrow=4, scales = 'free_y') +
  labs(x= 'Number of OTUs', y= 'Distance values')
ggsave('thesis/plots/norare_rare_relabund_normal_beatdiveristy.png', dpi=600)

# Rarefaction curve 
# Transform otutab into shape for vegan::rarecurve 
otutab_pre = shared %>% 
  pivot_wider(names_from = 'name', values_from = 'value') %>%
  filter(substr(Group, 1, 1) == 'M') %>%
  column_to_rownames('Group')
# Calculate rarefaction curve 
rarecurve(otutab_pre, step = 1000, xlab= 'Sample Size', ylab='OTUs')
ggsave('thesis/plots/rarefaction_curve.png', dpi=600)

# Rarefaction to 150 000 reads per sample with 999 iterattions. Which samples have less than 150000 redas?
insuff = rownames(otutab_pre)[rowSums(otutab_pre) < 150000]

# Exclude this samples from analysis

# Additional ideas
# Graph x= OTUs, y=mean relabund from most to least abundant
oturel = decostand(otutab_pre, MARGIN=1, method='total')
rowSums(oturel)

oturel %>%
  pivot_longer(names_to = 'name', values_to = 'value', cols= starts_with('Otu')) %>%
  group_by(name) %>%
  summarize(value = mean(value)) %>%
  ungroup() %>%
  ggplot(aes(x=name, y=value)) +
  geom_point(size = 1) +
  scale_y_log10() +
  labs(x = 'OTUs', y='log(Mean relative abundance)') +
  theme(axis.text.x=element_blank())
ggsave('thesis/plots/relabundanceOTUs.png', dpi=600)

##
# Prepare otutab, taxtab and metadata 
## 

# Min number of reads in sample and min sequences per OTU as determined in the exploration analysis
reads_per_sample = 150000 # we loose 6 samples
min_seqs_per_otu = sum(shared$value)*0.000001 # remove 0.0001% which is 0.000001

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

# Data for EtOH-EMA fraction VS microbiota samples 
otutabEM_pre = shared_pre %>%
  pivot_wider(names_from = 'name', values_from = 'value', values_fill = 0) %>%
  column_to_rownames('Group')

# Rarefy the data once - as we sequenced so deep that for the first analysis this is not crucial !
otutabEM = rrarefy(otutabEM_pre, sample=reads_per_sample)
saveRDS(otutabEM, 'thesis/data/otutab_ethanol_bulk.RDS')

# Extract OTUs that are present rarefied table 
otu_names = as.data.frame(otutabEM) %>% colnames() 

# Import taxonomy table
taxtab = read_tsv('data/mothur/final.opti_mcc.0.03.cons.taxonomy') %>%
  filter(OTU %in% otu_names) %>%
  select(name = "OTU", taxonomy = "Taxonomy") %>%
  mutate(taxonomy = str_replace_all(taxonomy, "\\\\|\\\"|\\(\\d+\\)", ""),
         taxonomy = str_replace(taxonomy, ";$", "")) %>%
  separate(taxonomy, into=c("Domain", "Phylum", "Class", "Order", "Family", "Genus"),
           sep=";") %>%
  mutate(Phylum = ifelse(Phylum %in% c('Firmicutes', 'Bacteroidetes', 'Actinobacteria', 'Proteobacteria', 'Bacteria_unclassified', 'Verrucomicrobia'), Phylum, 'Other')) %>%
  mutate(Phylum = case_when(
    Phylum == 'Firmicutes' ~ 'Bacillota',
    Phylum == 'Bacteroidetes' ~ 'Bacteroidota',
    Phylum == 'Actinobacteria' ~ 'Actinomycetota',
    Phylum == 'Proteobacteria' ~ 'Pseudomonadota',
    Phylum == 'Verrucomicrobia' ~ 'Verrucomicrobiota', 
    Phylum == 'Bacteria_unclassified' ~ 'unclassified Bacteria',
    TRUE ~ Phylum ))

saveRDS(taxtab,'thesis/data/taxtab.RDS')

# Import metadata
metadata = as_tibble(read.csv('data/metadata.csv', sep=';')) %>%
  mutate(date=dmy(date)) %>%
  filter(Group %in% rownames(otutabEM)) %>%
  mutate(biota = ifelse(biota == 'microbiota', 'Bulk microbiota', 'Ethanol treated sample'))
saveRDS(metadata, 'thesis/data/metadata.RDS')

## Normalized abundances from ddPCR experiment 
# function to read in Quantasoft results
read.quantasoft = function(file) {
  # Read the lines
  lines = readLines(file)
  # Extract header 
  header = strsplit(lines[1], ',')[[1]]
  # Spilt by , 
  data_lines = sapply(lines[-1], function(line) strsplit(substr(line, 2, nchar(line)), ",")[[1]])
  # turn into data.frame and add header and remove row.names
  data_df = as.data.frame(t(data_lines), stringsAsFactors = FALSE)
  names(data_df) = header
  rownames(data_df) = NULL
  # Select only the info I need, remove " from names of well and samples
  data_df_fin = data_df %>% select(Well, Sample, Concentration, CopiesPer20uLWell, Positives, Negatives, AcceptedDroplets) %>%
    filter(Concentration != 'No Call') %>%
    transform(Concentration = as.numeric(Concentration), 
              CopiesPer20uLWell = as.numeric(CopiesPer20uLWell), 
              Positives = as.numeric(Positives),
              Negatives = as.numeric(Negatives),
              AcceptedDroplets = as.numeric(AcceptedDroplets)) %>%
    mutate(across(c(Well, Sample), ~gsub('"', '', .)))
  
}

# 
samples_info = read_excel('data/vzorci.xlsx', sheet = 6)

plate1_m  = read.quantasoft('data/absolute_quantification/original_files/20240322_ddPCR_v3v4_microbiota_plate1_results.csv') %>%
  filter(AcceptedDroplets > 10000 & Positives > 12)

plate2_m = read.quantasoft('data/absolute_quantification/original_files/20240322_ddPCR_v3v4_microbiota_plate2_results.csv') %>%
  filter(AcceptedDroplets > 10000 & Positives > 8)

plates_e = read.quantasoft('data/absolute_quantification/original_files/20240513_ddPCR_v3v4_sporobiota_1_results.csv') %>%
  rbind(read.quantasoft('data/absolute_quantification/original_files/20240513_ddPCR_v3v4_sporobiota_2_results.csv')) %>%
  filter(AcceptedDroplets > 10000 & Positives > 1) %>%
  # Exclude from analysis because the amount of DNA was insufficient SB008, SB009, SE002, SE003, SF001, SF009, SH007, SC013
  filter(!(Sample %in% c('SB008', 'SB009', 'SE002', 'SE003', 'SF001', 'SF009', 'SH007', 'SC013')))

sample = read.quantasoft('data/absolute_quantification/original_files/MA001.csv') 

ddPCR = rbind(plate1_m, plate2_m, plates_e, sample) %>%
  group_by(Sample) %>%
  summarise(Concentration = mean(Concentration)) %>%
  left_join(samples_info, by =join_by('Sample'=='Group')) %>%
  # Calculate the copy number of 16s rRNA gene per ng of DNA
  # concentration = Poisson correlted value copies per ul
  # 25/2.5 = adjust for the amount of DNA in reaction
  # 25/20 = adjust for the reaction made VS reaction used
  # dilution of the DNA 
  # dilution from original samples to normalized value
  mutate(copies = (Concentration * (25/2.5) * (25/20) * dilution_ddPCR * (DNAconc/DNAconc_seq)))
saveRDS(ddPCR, 'thesis/data/ddPCR.RDS')

# Calculate normalized abundances 
# Multiply relative abundances by CopiesPerngDNA = absolute abundance per ng DNA OR 
# CopiesPerulSample = absolutne abundance per ul DNA in original sample

otutab_normrel = rownames_to_column(as.data.frame(otutabEM), 'Group') %>%
  pivot_longer(cols = starts_with('Otu')) %>%
  group_by(Group) %>%
  mutate(rel_abund = value/sum(value)) %>%
  ungroup() %>%
  left_join(ddPCR, by = join_by('Group' == 'Sample')) %>%
  mutate(norm_abund = rel_abund*copies) %>%
  select(Group, name, value, rel_abund, norm_abund) %>%
  filter(!is.na(norm_abund))

saveRDS( otutab_normrel, 'thesis/data/otutab_normrel.RDS')
