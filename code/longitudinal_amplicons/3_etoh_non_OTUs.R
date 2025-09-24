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

otutab <- readRDS('longitudinal_amplicons/data/otutab_ethanol_bulk.RDS')
taxtab <- readRDS('longitudinal_amplicons/data/taxtab.RDS')
metadata <- readRDS('longitudinal_amplicons/data/metadata.RDS')
norm_rel <- readRDS('longitudinal_amplicons/data/otutab_normrel.RDS') %>%
  left_join(select(metadata, Group, original_sample), by = 'Group')

#
# Colors to be used
col2 <- c('#66a5f1', '#f18e66')
col4 <- c('#f0a336', '#3CB371', '#f35020', '#6490f0')

# Define ethanol resistant OTUs and seqs! 
otu_long <- norm_rel %>% 
  mutate(PA = ifelse(value > 0, 1, 0)) %>%
  left_join(metadata %>% select(Group, person, date), by = 'Group') %>%
  left_join(taxtab, by = 'name')

# OTU that is ethanol resistant once is always ethanol resistant 
etoh_otus <- left_join(otu_long %>% filter(substr(Group, 1, 1) == 'M'), 
                       otu_long %>% filter(substr(Group, 1, 1) == 'S'), 
                       by = join_by('name', 'original_sample', 'person', 'Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus')) %>%
  # Define if OTU in a sample of stool is ethanol resistant 
  # Contition 1: present in both bulk microbiota sample and ethanol resistant fraction
  # Condition 2: higher relative abudnance in EtOH sample than microbiota
  mutate(is_etoh_resistant = ifelse(value.x > 0 & value.y > 0 & rel_abund.y > rel_abund.x, 'Yes', 'No')) %>%
  group_by(name) %>%
  # Calculate the number of times this OTU was present in samples
  reframe(no_present = sum(PA.x), 
          # Caluclate how many times OTU was defined as part of EtOH fraction based on Conditions 1 & 2
          no_Yes = ceiling(sum(is_etoh_resistant == 'Yes', na.rm = TRUE))) %>%
  ungroup() %>%
  # OTUs that have been defined as part Of the ethanol resistant fraction in at least 5% of samples where they were found! 
  # (to avoid mistakes of protocol and exclude highly abundant OTUs that maybe were seen as ethanol resistant but just didn't get destoryed!)
  filter(no_Yes > (no_present * 0.05)) %>%
  #filter(no_Yes > (no_present * 0.1)) %>%
  #filter(no_Yes > 1) %>%
  # Extract names! 
  pull(unique(name))

length(unique(etoh_otus))

uncertain_otus <- left_join(otu_long %>% filter(substr(Group, 1, 1) == 'M'), 
                            otu_long %>% filter(substr(Group, 1, 1) == 'S'), 
                            by = join_by('name', 'original_sample', 'person', 'Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus')) %>%
  mutate(is_etoh_resistant = ifelse(value.x > 0 & value.y > 0 & rel_abund.y > rel_abund.x, 'Yes', 'No')) %>%
  group_by(name) %>%
  reframe(no_present = sum(PA.x), 
          no_Yes = ceiling(sum(is_etoh_resistant == 'Yes', na.rm = TRUE))) %>%
  ungroup() %>%
  # Filter OTUs that were detected as EtOH resistant at least once, but were detected as such in less than 5% of samples, to exclude them from the analysis 
  filter(no_Yes > 1) %>%
  filter(no_Yes < (no_present * 0.05)) %>%
  pull(unique(name))
length(unique(uncertain_otus))

nonetoh_otus <- otu_long %>% filter(substr(Group, 1, 1) == 'M' & PA == 1) %>%
  filter(!(name %in% uncertain_otus) & !(name %in% etoh_otus)) %>%
  pull(unique(name))
length(unique(nonetoh_otus))


##  OTU analysis 

# Create the 4 fractions 
# Ethanol resistant OTUs AND non-ethanol resistant OTUs + divide by phylum (Bacillota + other)
# At the level of Bacillota 
etoh_bacillota <- filter(otu_long, substr(Group, 1, 1) == 'M' & name %in% etoh_otus & Phylum == 'Bacillota') %>%
  mutate(Group = paste0(Group, "-EB"), is_ethanol_resistant = 'Ethanol resistant', taxonomy = 'Bacillota', fraction = 'Ethanol resistant Bacillota')
# min = 78

non_etoh_bacillota <-  filter(otu_long, substr(Group, 1, 1) == 'M' & name %in% nonetoh_otus & Phylum == 'Bacillota') %>%
  mutate(Group = paste0(Group, "-NB"), is_ethanol_resistant = 'Ethanol non-resistant', taxonomy = 'Bacillota', fraction = 'Ethanol non-resistant Bacillota')
# min = 343

etoh_other <- filter(otu_long, substr(Group, 1, 1) == 'M' & name %in% etoh_otus & Phylum != 'Bacillota') %>%
  mutate(Group = paste0(Group, "-E"), is_ethanol_resistant = 'Ethanol resistant', taxonomy = 'Other taxa', fraction = 'Other ethanol resistant taxa') 
# min = 24

non_etoh_other <- filter(otu_long, substr(Group, 1, 1) == 'M' & name %in% nonetoh_otus & Phylum != 'Bacillota') %>% 
  mutate(Group = paste0(Group, "-NE"), is_ethanol_resistant = 'Ethanol non-resistant', taxonomy = 'Other taxa', fraction = 'Other ethanol non-resistant taxa')
# min = 64

##
long_all <- rbind(etoh_bacillota, non_etoh_bacillota, etoh_other, non_etoh_other)

long_all$Phylum <- factor(long_all$Phylum, levels = c('Bacillota', 'Bacteroidota', 'Actinomycetota', 'Pseudomonadota', 
                                                      'Verrucomicrobiota', 'unclassified Bacteria', 'Other'))
saveRDS(long_all, 'longitudinal_amplicons/data/long_all.RDS')

# The number and relative abundance of ethanol resistant fraction within bulk microbiota samples 
res_relative <- data.frame()
for (i in unique(long_all$Phylum)) {
  sub <- filter(long_all, Phylum == i)
  res <- kruskal.test(sub$rel_abund, sub$is_ethanol_resistant)
  res_relative <- rbind(res_relative, data.frame(Phylum = i, 
                                                 pvalue = res$p.value, 
                                                 statistic = res$statistic))
}

res_relative

relative <- ggplot(long_all) +
  geom_boxplot(mapping = aes(x = Phylum, y = rel_abund, fill = is_ethanol_resistant)) +
  geom_text(res_relative, mapping =aes(x = Phylum, y = 1, label = paste('p =', scientific(pvalue, digits =0))), size = 4) +
  scale_y_log10() +
  scale_fill_manual(values = col2) +
  labs(x = '', y = 'log10(relative abundance)', fill = '') +
  theme(legend.position = 'bottom', 
        axis.ticks.x = element_blank(), 
        plot.margin = unit(c(0, 0.2, 0.2, 0.2), "cm")) +
  guides(fill = guide_legend(ncol = 4))
relative

# The number of unique OTUs in each phylum 
number <- long_all %>%
  group_by(is_ethanol_resistant, Phylum) %>%
  reframe(no_otus = n_distinct(name), 
          sum_value = sum(value)) %>%
  group_by(is_ethanol_resistant) %>%
  mutate(per = sum_value / sum(sum_value)* 100) %>%
  ggplot(aes(x = Phylum, y = no_otus, fill = is_ethanol_resistant)) +
  geom_col(position = position_dodge()) +
  geom_text(aes(label = paste(no_otus),
                vjust = ifelse(no_otus > 1000, 1.1, - 0.3)), size = 4,
            position = position_dodge(width = 0.9)) +
  scale_fill_manual(values = col2) +
  # coord_cartesian(ylim = c(0, 1700)) +
  labs(x = '', y = 'Number of OTUs', fill = '') +
  theme(legend.position = 'bottom', 
        # margin(t, r, l, b)
        plot.margin = unit(c(0.1, 0.2, 0.2, 0.1), "cm")) + 
  guides(fill = guide_legend(ncol = 4))

number
ggarrange(relative + labs(tag = 'A'), number + labs(tag = 'B'), 
          common.legend = TRUE, legend = 'bottom', ncol = 1, align = 'v', heights = c(0.8, 1))
ggsave('longitudinal_amplicons/plots/realtive_number_etoh_non_OTUs.png', dpi = 600)

# 
# Are ethanol resistant OTUs more likely to be shared or present in a single individual? 
# An OTU is present in an individual, if we saw it in at elast 1/3 of the samples (n=4).

# Beta diversity 

# Minimal number of OTUs in a fraction 
calculate_min <- function(otu_data) {
  min <-  otu_data %>%
    group_by(Group) %>%
    summarise(sum = sum(PA), .groups = 'drop') %>%
    summarise(min = min(sum)) %>%
    pull(min)
}


etoh_bac_min <- calculate_min(etoh_bacillota) #78
etoh_other_min <- calculate_min(etoh_other) # 24
non_etoh_bac_min <- calculate_min(non_etoh_bacillota) # 292
non_etoh_other_min <- calculate_min(non_etoh_other) # 60

# Function to calculate beta distances (Bray-Curtis OR Jaccard)
min <- 24
calculate_dist <- function(otu_data, method) {
  dist_all <- data.frame()
  
  meta <- distinct(otu_data, Group, person, date, fraction, is_ethanol_resistant, taxonomy)
  
  # min <- otu_data %>%
  #   group_by(Group) %>%
  #   summarise(sum = sum(PA), .groups = 'drop') %>%
  #   summarise(min = min(sum)-5) %>%
  #   pull(min)
  
  otutab <- otu_data %>%
    select(Group, name, value) %>%
    pivot_wider(names_from = 'name', values_from = 'value', values_fill = 0) %>%
    column_to_rownames('Group')
  
  for (i in 1:999) {
    # Resample OTUs within each fraction
    otutab_t <- t(otutab)
    resampled_t <- otutab_t[sample(1:nrow(otutab_t), size = min, replace = TRUE), ]
    resampled_otutab <- t(resampled_t)
    
    # Calculate distances (Bray-Curtis)
    dist <- vegdist(resampled_otutab, method = method)
    
    # Tidy the Bray-Curtis matrix
    dist_long <- as.matrix(dist) %>%
      as_tibble(rownames = 'Group') %>%
      pivot_longer(-Group) %>%
      filter(Group != name)
    
    dist_all <- rbind(dist_all, dist_long)
  }
  
  dist <- dist_all %>%
    mutate(sample_pairs = paste(Group, name)) %>%
    group_by(sample_pairs) %>%
    summarise(mean_value = mean(value, na.rm = TRUE), 
              median_value = median(value, na.rm = TRUE),
              sd = sd(value, na.rm = TRUE), .groups = 'drop') %>%
    ungroup() %>%
    separate(sample_pairs, into = c("Group", "name"), sep = " ") %>%
    left_join(meta, by = 'Group') %>%
    left_join(meta, by = join_by('name' == 'Group', 'fraction', 'is_ethanol_resistant', 'taxonomy')) %>%
    mutate(same_person = ifelse(person.x == person.y, 'Intra individual', 'Inter individual'), 
           date_dist = abs(date.x-date.y))
  
  return(dist)
}

# Calculate Bray-Curtis distances and combine all results 
bray <- calculate_dist(etoh_bacillota, 'bray') %>%
  rbind(calculate_dist(non_etoh_bacillota, 'bray')) %>%
  rbind(calculate_dist(etoh_other, 'bray')) %>%
  rbind(calculate_dist(non_etoh_other, 'bray'))

# Kruskal test for Bray-Curtis distances 
# Within individual 
bray_within <- filter(bray, same_person == 'Intra individual')

kruskal_within <- kruskal.test(median_value~fraction, data = bray_within)
# But between which groups ?
wilcox_within <- pairwise.wilcox.test(bray_within$median_value, bray_within$fraction, paired = FALSE)
# Between 
bray_between <- filter(bray, same_person == 'Inter individual')
kruskal_between <- kruskal.test(median_value ~fraction, data = bray_between)
wilcox_between <- pairwise.wilcox.test(bray_between$median_value, bray_between$fraction, paired = FALSE)

wilcox_to_df <- function(wilcox_result, same_person_label) {
  # Extract the matrix of p-values
  pval <- as.data.frame(as.table(wilcox_result$p.value)) %>%
    na.omit()
  
  # Rename columns for clarity
  colnames(pval) <- c("fraction1", "fraction2", "pvalue")
  
  # Add the same_person column
  pval$same_person <- same_person_label
  
  return(pval)
}

wilcox_bray <- rbind(wilcox_to_df(wilcox_between, 'Inter individual'), 
                     wilcox_to_df(wilcox_within, 'Intra individual')) %>%
  mutate(is_ethanol_resistant = ifelse(fraction1 %in% c('Ethanol resistant Bacillota', 'Other ethanol resistant taxa'), 'Ethanol resistant', 'Ethanol non-resistant'), 
         taxonomy = ifelse(fraction1 %in% c('Ethanol resistant Bacillota', 'Ethanol non-resistant Bacillota'), 'Bacillota', 'Other taxa'), 
         is_ethanol_resistant2 = ifelse(fraction2 %in% c('Ethanol resistant Bacillota', 'Other ethanol resistant taxa'), 'Ethanol resistant', 'Ethanol non-resistant'), 
         taxonomy2 = ifelse(fraction2 %in% c('Ethanol resistant Bacillota', 'Ethanol non-resistant Bacillota'), 'Bacillota', 'Other taxa')) %>%
  filter(taxonomy == taxonomy2) %>%
  select(fraction1, fraction2, pvalue, same_person, is_ethanol_resistant, taxonomy)

# Correlations for distance / time 
time_corr <- function(data) {
  # Function to compute correlation and p-value for each subset
  cor_function <- function(df) {
    cor_result <- cor.test(as.numeric(df$date_dist), df$median, method = "pearson")
    return(data.frame(corr = cor_result$estimate, p_value = cor_result$p.value))
  }
  
  # Apply the cor_function for each unique fraction and store the results in a data frame
  corr_table <- data %>%
    group_by(fraction) %>%
    do({
      cor_function(.)
    }) %>%
    ungroup()
  
  return(corr_table)
}

time_bray <- bray %>%
  # Filter different individuals
  filter(same_person == 'Intra individual') %>%
  # group by difference between days and person
  group_by(fraction, is_ethanol_resistant, taxonomy,  date_dist) %>%
  reframe(median=median(median_value), sd= sd(median_value)) %>%
  ungroup()


# Calculate correaltions between diff (time between samplings) and distance metric
bray_corr_time <- time_corr(time_bray)
bray_corr_time

bray_boxplot <- ggplot(bray) +
  geom_boxplot(mapping = aes(x = taxonomy, y = median_value, fill = is_ethanol_resistant)) +
  geom_line(mapping = aes(x = .25, y = .25, linetype = taxonomy)) +
  geom_text(data = wilcox_bray, mapping = aes(y = .05, x = taxonomy, label = paste('p =', scientific(pvalue, digits = 0)))) + 
  scale_fill_manual(values = col) +
  labs(y = 'Bray-Curtis dissimilarity', x = '', fill = '', linetype = 'Phylum') +
  theme(legend.position = 'bottom', axis.ticks.x = element_blank()) +
  guides(fill = guide_legend(ncol = 4)) +
  facet_grid(~same_person) 

b_time <- time_bray %>%
  ggplot(aes(x = date_dist, y = median, color = is_ethanol_resistant)) +
  geom_point() +
  geom_smooth(method = 'lm', se = TRUE, alpha = .2, mapping = aes(linetype = taxonomy, color = is_ethanol_resistant)) +
  scale_color_manual(values = col) +
  labs(x = 'Days between sampling', y = 'Bray-Curtis dissimilarity', color = '', linetype = 'Phylum') +
  theme(legend.position = 'bottom')
b_time

##
# Jaccard 
jaccard <- calculate_dist(etoh_bacillota, 'jaccard') %>%
  rbind(calculate_dist(non_etoh_bacillota, 'jaccard')) %>%
  rbind(calculate_dist(etoh_other, 'jaccard')) %>%
  rbind(calculate_dist(non_etoh_other, 'jaccard'))

# Kruskal test for jaccard distances 
# Within individual 
jaccard_within <- filter(jaccard, same_person == 'Intra individual')

kruskal_within <- kruskal.test(median_value~fraction, data = jaccard_within)
# But between which groups ?
wilcox_within <- pairwise.wilcox.test(jaccard_within$median_value, jaccard_within$fraction, paired = FALSE)
# Between 
jaccard_between <- filter(jaccard, same_person == 'Inter individual')
kruskal_between <- kruskal.test(median_value ~fraction, data = jaccard_between)
wilcox_between <- pairwise.wilcox.test(jaccard_between$median_value, jaccard_between$fraction, paired = FALSE)

wilcox_jaccard <- rbind(wilcox_to_df(wilcox_between, 'Inter individual'), 
                        wilcox_to_df(wilcox_within, 'Intra individual')) %>%
  mutate(is_ethanol_resistant = ifelse(fraction1 %in% c('Ethanol resistant Bacillota', 'Other ethanol resistant taxa'), 'Ethanol resistant', 'Ethanol non-resistant'), 
         taxonomy = ifelse(fraction1 %in% c('Ethanol resistant Bacillota', 'Ethanol non-resistant Bacillota'), 'Bacillota', 'Other taxa'), 
         is_ethanol_resistant2 = ifelse(fraction2 %in% c('Ethanol resistant Bacillota', 'Other ethanol resistant taxa'), 'Ethanol resistant', 'Ethanol non-resistant'), 
         taxonomy2 = ifelse(fraction2 %in% c('Ethanol resistant Bacillota', 'Ethanol non-resistant Bacillota'), 'Bacillota', 'Other taxa')) %>%
  filter(taxonomy == taxonomy2) %>%
  select(fraction1, fraction2, pvalue, same_person, is_ethanol_resistant, taxonomy)

# Correlations for distance / time 
time_jaccard <- jaccard %>%
  # Filter different individuals
  filter(same_person == 'Intra individual') %>%
  # group by difference between days and person
  group_by(fraction, is_ethanol_resistant, taxonomy,  date_dist) %>%
  reframe(median=median(median_value), sd= sd(median_value)) %>%
  ungroup()


# Calculate correaltions between diff (time between samplings) and distance metric
jaccard_corr_time <- time_corr(time_jaccard)
jaccard_corr_time

jaccard_boxplot <- ggplot(jaccard) +
  geom_boxplot(mapping = aes(x = taxonomy, y = median_value, fill = is_ethanol_resistant)) +
  geom_line(mapping = aes(x = .25, y = .25, linetype = taxonomy)) +
  geom_text(data = wilcox_jaccard, mapping = aes(y = .05, x = taxonomy, label = paste('p =', scientific(pvalue, digits = 0)))) + 
  scale_fill_manual(values = col) +
  labs(y = 'Jaccard distance', x = '', fill = '', linetype = 'Phylum') +
  theme(legend.position = 'bottom', axis.ticks.x = element_blank()) +
  guides(fill = guide_legend(ncol = 4)) +
  facet_grid(~same_person) 

j_time <- time_jaccard %>%
  ggplot(aes(x = date_dist, y = median, color = is_ethanol_resistant)) +
  geom_point() +
  geom_smooth(method = 'lm', se = TRUE, alpha = .2, mapping = aes(linetype = taxonomy, color = is_ethanol_resistant)) +
  scale_color_manual(values = col) +
  labs(x = 'Days between sampling', y = 'Jaccard distance', color = '', linetype = 'Phylum') +
  theme(legend.position = 'bottom')

# Combine plots with a shared legend
ggarrange(bray_boxplot + labs(tag = 'A'), 
          b_time + labs(tag = 'B'), common.legend = TRUE, legend = 'bottom',ncol=2, widths = c(0.8, 1))

ggsave('endospore_dynamics/out/supplement_figure1.png', dpi = 600)

ggarrange(jaccard_boxplot + labs(tag = 'A'), j_time + labs(tag = 'B'),
          nrow = 1, ncol = 2, common.legend = TRUE, legend = 'bottom')
ggsave('endospore_dynamics/out/figure1.png', dpi=600)


ggarrange(bray_boxplot + labs(tag = 'A'), b_time + labs(tag = 'B'), 
          ncol = 2, common.legend = TRUE, legend = 'bottom')
ggsave('endospore_dynamics/out/supplement_figure2.png', dpi=600)


# Additional test for usage of beta diveristy metrics: Is the difference we see between ethanol resistant and ethanol non-resistant only, 
# becouse of different relative abundances of OTUs (Bacillota have higher relative abundance)

sf <- long_all %>%
  group_by(name) %>%
  reframe(mean_rel_abund =  mean(rel_abund), 
          sumsq_diff_abund = sum((outer(rel_abund, rel_abund, `-`)^2)[lower.tri(outer(rel_abund, rel_abund))])) %>%
  left_join(select(long_all, name, fraction, is_ethanol_resistant), by = 'name')

ggplot(sf, aes(x = log10(mean_rel_abund), y = log10(sumsq_diff_abund), color = is_ethanol_resistant)) +
  geom_point() +
  scale_color_manual(values = col) +
  labs(x = 'Mean relative abundance of OTU', y = 'Sum of squared differences between realtive abudnances of an OTU in different samples', color = '')
ggsave('endospore_dynamics/out/supplement_figure8.png', dpi = 600)


ggplot(sf, aes(x = log10(mean_rel_abund), y = log10(sumsq_diff_abund), color = fraction)) +
  geom_point() +
  scale_color_manual(values = col4) +
  labs(x = 'Mean relative abundance of OTU', y = 'Sum of squared differences between realtive abudnances of an OTU in different samples', color = '') 
ggsave('endospore_dynamics/out/supplement_figure8_fractions1.png', dpi = 600)


