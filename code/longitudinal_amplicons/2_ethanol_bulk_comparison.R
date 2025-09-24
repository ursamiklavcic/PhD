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

col_phylum = c('#1F77B4', '#FF7F0E',  '#2CA02C',  '#D62728', '#9467BD', '#8C564B', '#f4d03f')
cole = c('#d87328')
colm = c('#27ae60')
colme = c('#27ae60', '#d87328' )

# efficiency of ethanol treatment on stool samples in this study
norm_rel2 <- filter(norm_rel, substr(Group, 1, 1) == 'M') %>%
  left_join(filter(norm_rel, substr(Group, 1, 1) == 'S'), by = c('name', 'original_sample')) %>%
  left_join(taxtab, by = 'name')

norm_rel2 %>%
  filter(rel_abund.x != 0 & rel_abund.y != 0) %>%
  mutate(ratio = rel_abund.x/rel_abund.y) %>%
  group_by(original_sample, Phylum) %>%
  reframe(mean_ratio = mean(ratio)) %>%
  ggplot(aes(x = mean_ratio, y = Phylum)) +
  geom_boxplot() +
  scale_x_log10() +
  geom_point(size = 2, alpha = .2) +
  geom_vline(xintercept = 1) +
  labs(x = 'Ratio of relative abundance in stool samples vs ethanol treated sample', y = '') 
ggsave('longitudinal_amplicons/plots/EMenrichment_phylum.png', dpi=600)

# Data frame of everthing! 
otu_all <- norm_rel %>%
  left_join(metadata, by = 'Group') %>%
  left_join(taxtab, by = 'name')

otu_all$Phylum <- factor(otu_all$Phylum, levels = c('Bacillota', 'Bacteroidota', 'Actinomycetota', 'Pseudomonadota',
                                                    'Verrucomicrobiota', 'unclassified Bacteria', 'Other'))
# Relative abundance plots for comparison
otu_all %>%
  group_by(biota) %>% 
  mutate(rel = rel_abund/sum(rel_abund)) %>%
  ungroup() %>%
  ggplot(aes(x = biota, y = rel, fill = Phylum)) +
  geom_col() +
  scale_fill_manual(values = col_phylum) +
  scale_y_continuous(labels = function(y) paste0(y*100, "%")) +
  labs(x = '', y= 'Relative abundance', fill = 'Phylum')
ggsave('longitudinal_amplicons/plots/EM_rel_abund.png', dpi=600)

# Normalized abundance plots for comparison 
ggarrange(otu_all %>% 
            filter(biota == 'Bulk microbiota') %>%
            ggplot(aes(x = biota, y = norm_abund, fill = Phylum)) +
            geom_col() +
            scale_fill_manual(values = col_phylum) +
            labs(x = '', y= 'Normalized abundance', fill = ''), 
          otu_all %>% 
            filter(biota == 'Ethanol treated sample') %>%
            ggplot(aes(x = biota, y = norm_abund, fill = Phylum)) +
            geom_col() +
            scale_fill_manual(values = col_phylum) +
            labs(x = '', y= 'Normalized abundance', fill = ''), 
          common.legend = TRUE, legend = 'bottom')

ggsave('longitudinal_amplicons/plots/EM_norm_abund.png', dpi=600)

# 
## Alpha diveristy 
# 
richness <- estimateR(otutab) # observed richness and Chao1
evenness <- diversity(otutab)/log(specnumber(otutab)) # evenness index
shannon <- diversity(otutab, index = 'shannon')

# Join all calculations and metadata
alpha_meta <- as_tibble(as.list(evenness)) %>% pivot_longer(names_to = 'Group', values_to = 'evenness', cols = starts_with(c('M', 'S'))) %>%
  left_join(t(richness) %>% as.data.frame() %>% rownames_to_column('Group'), by='Group') %>%
  left_join(as_tibble(as.list(shannon)) %>% pivot_longer(names_to = 'Group', values_to = 'shannon', cols = starts_with(c('M', 'S')))) %>%
  left_join(metadata, by='Group') %>%
  mutate(person2 = person)

# Observed number of OTUs
# between hosts
ggplot(alpha_meta, aes(x = person, y = S.obs, color = person )) +
  geom_boxplot() +
  geom_jitter(alpha=.2, width = .1) +
  facet_wrap(~ biota, nrow = 1) +
  labs(x = 'Individuals', y = '#OTUs') +
  theme(legend.position = 'none')
ggsave('longitudinal_amplicons/plots/EM_observed_between_hosts.png', dpi=600)

# in time within host 
ggplot(alpha_meta, aes(x=day, y=S.obs)) +
  geom_line(data=alpha_meta %>% dplyr::select(-person) %>% filter(biota == 'Bulk microbiota'), 
            aes(group=person2), color= colm, linewidth=0.5, alpha=0.5) +
  geom_line(data=alpha_meta %>% dplyr::select(-person) %>% filter(biota == 'Ethanol treated sample'), 
            aes(group=person2), color= cole, linewidth=0.5, alpha=0.5)+
  geom_line(data=alpha_meta %>% filter(biota == 'Bulk microbiota'),
            aes(color=person), color= colm, linewidth=1.2) +
  geom_line(data=alpha_meta %>% filter(biota == 'Ethanol treated sample'), 
            color=cole, linewidth=1.2) +
  facet_wrap(~person, scales = 'free') +
  labs(x='Day', y= '#OTUs', color = 'Sample')
ggsave('longitudinal_amplicons/plots/EM_observed_time.png', dpi=600)

# Evenness
# between hosts
ggplot(alpha_meta, aes(x = person, y = evenness, color = person )) +
  geom_boxplot() +
  geom_jitter(alpha=.2, width = .1) +
  facet_wrap(~ biota, nrow = 1) +
  labs(x = 'Individuals', y = 'Evenness') +
  theme(legend.position = 'none')
ggsave('longitudinal_amplicons/plots/EM_evenness_between_hosts.png', dpi=600)

# in time within host 
ggplot(alpha_meta, aes(x=day, y=evenness)) +
  geom_line(data=alpha_meta %>% dplyr::select(-person) %>% filter(biota == 'Bulk microbiota'), 
            aes(group=person2), color= colm, linewidth=0.5, alpha=0.5) +
  geom_line(data=alpha_meta %>% dplyr::select(-person) %>% filter(biota == 'Ethanol treated sample'), 
            aes(group=person2), color= cole, linewidth=0.5, alpha=0.5)+
  geom_line(data=alpha_meta %>% filter(biota == 'Bulk microbiota'),
            aes(color=person), color= colm, linewidth=1.2) +
  geom_line(data=alpha_meta %>% filter(biota == 'Ethanol treated sample'), 
            color=cole, linewidth=1.2) +
  facet_wrap(~person, scales = 'free') +
  labs(x='Day', y= 'Evenness', color = 'Sample')
ggsave('longitudinal_amplicons/plots/EM_evenness_time.png', dpi=600)

# Shannon
# between hosts
ggplot(alpha_meta, aes(x = person, y = shannon, color = person )) +
  geom_boxplot() +
  geom_jitter(alpha=.2, width = .1) +
  facet_wrap(~ biota, nrow = 1) +
  labs(x = 'Individuals', y = 'Shannon') +
  theme(legend.position = 'none')
ggsave('longitudinal_amplicons/plots/EM_shannon_between_hosts.png', dpi=600)

# in time within host 
ggplot(alpha_meta, aes(x=day, y=shannon)) +
  geom_line(data=alpha_meta %>% dplyr::select(-person) %>% filter(biota == 'Bulk microbiota'), 
            aes(group=person2), color= colm, linewidth=0.5, alpha=0.5) +
  geom_line(data=alpha_meta %>% dplyr::select(-person) %>% filter(biota == 'Ethanol treated sample'), 
            aes(group=person2), color= cole, linewidth=0.5, alpha=0.5)+
  geom_line(data=alpha_meta %>% filter(biota == 'Bulk microbiota'),
            aes(color=person), color= colm, linewidth=1.2) +
  geom_line(data=alpha_meta %>% filter(biota == 'Ethanol treated sample'), 
            color=cole, linewidth=1.2) +
  facet_wrap(~person, scales = 'free') +
  labs(x='Day', y= 'Shannon', color = 'Sample')
ggsave('longitudinal_amplicons/plots/EM_shannon_time.png', dpi=600)

# Correlations between bulk microbiota and ethanol treated samples 
# Evenness correlation 
evenness <- alpha_meta %>% select(original_sample, biota, evenness, person) %>%
  pivot_wider(names_from = 'biota', values_from = 'evenness') %>%
  ggplot(aes(x = `Bulk microbiota`, y= `Ethanol treated sample`)) +
  geom_point() +
  geom_smooth(method = 'lm') +
  ggpubr::stat_cor(method = 'pearson', aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~"))) +
  labs(subtitle = 'Evenness')

observed <- alpha_meta %>% select(original_sample, biota, S.obs, person) %>%
  pivot_wider(names_from = 'biota', values_from = 'S.obs') %>%
  ggplot(aes(x = `Bulk microbiota`, y= `Ethanol treated sample`)) +
  geom_point() +
  geom_smooth(method = 'lm') +
  ggpubr::stat_cor(method = 'pearson', aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~"))) +
  labs(subtitle = '#OTUs')

shannon <- alpha_meta %>% select(original_sample, biota, shannon, person) %>%
  pivot_wider(names_from = 'biota', values_from = 'shannon') %>%
  ggplot(aes(x = `Bulk microbiota`, y= `Ethanol treated sample`)) +
  geom_point() +
  geom_smooth(method = 'lm') +
  ggpubr::stat_cor(method = 'pearson', aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~"))) +
  labs(subtitle = 'Shannon')

ggarrange(evenness, observed, shannon, 
          nrow = 1)
ggsave('longitudinal_amplicons/plots/EM_alpha_correlation.png', dpi=600)


# Beta diversity 
# As samples were isoalted with 2 different protocols, we should not calculate beta diversity for all of them, but each protocol sperately 
calc_dist <- function(tab, method, i ) {
  # Separate the smaples 
  tab_filt <- as.data.frame(tab) %>%
    rownames_to_column('Group') %>%
    filter(substr(Group, 1, 1) == i) %>%
    column_to_rownames('Group') %>%
    as.matrix()
    
  # Calculate dissimilarity / distance
  dist <- vegdist(tab_filt, method = method)
  
  # Tidy the matrix 
  dist_long <- as.matrix(dist) %>%
    as_tibble(rownames = 'Group') %>%
    pivot_longer(-Group) %>%
    filter(Group != name) %>%
    # Add metadata 
    left_join(metadata %>% select(Group, person, date), by = 'Group') %>%
    left_join(metadata %>% select(Group, person, date), by = join_by('name' == 'Group')) %>%
    # Define inter/intra inidivdual etc. 
    mutate(same_person = ifelse(person.x == person.y, 'Intra individual', 'Inter individuals'), 
           date_dist = abs(date.x - date.y))
  
  return(dist_long)
}

bray <- calc_dist(otutab, 'bray', 'M') %>%
  mutate(biota = 'Bulk microbiota') %>%
  rbind(calc_dist(otutab, 'bray', 'S') %>% 
          mutate(biota = 'Ethanol treated sample'))

# Kruskal test for Bray-Curtis distances 
# Within individual 
bray_within <- filter(bray, same_person == 'Intra individual')
kruskal_within <- kruskal.test(value~biota, data = bray_within)
# But between which groups ?
wilcox_within <- pairwise.wilcox.test(bray_within$value, bray_within$biota, paired = FALSE)

# Between 
bray_between <- filter(bray, same_person == 'Inter individuals')
kruskal_between <- kruskal.test(value ~ biota, data = bray_between)
wilcox_between <- pairwise.wilcox.test(bray_between$value, bray_between$biota, paired = FALSE)

wilcox_bray <- data.frame(pvalue = wilcox_within$p.value, same_person = 'Intra individual') %>%
  rbind(data.frame(pvalue = wilcox_between$p.value, same_person = 'Inter individuals'))
colnames(wilcox_bray) <- c('pvalue', 'same_person')

# Correlations for distance / time 
time_corr <- function(data) {
  # Function to compute correlation and p-value for each subset
  cor_function <- function(df) {
    cor_result <- cor.test(as.numeric(df$date_dist), df$median, method = "pearson")
    return(data.frame(corr = cor_result$estimate, p_value = cor_result$p.value))
  }
  
  # Apply the cor_function for each unique fraction and store the results in a data frame
  corr_table <- data %>%
    group_by(biota) %>%
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
  group_by(biota, date_dist) %>%
  reframe(median=median(value)) %>%
  ungroup()

# Calculate correaltions between diff (time between samplings) and distance metric
bray_corr_time <- time_corr(time_bray)
bray_corr_time

bray_boxplot <- ggplot(bray) +
  geom_boxplot(mapping = aes(x = same_person, y = value, fill = biota)) +
  geom_text(wilcox_bray, mapping = aes( x= same_person, y = 1, label = paste('p=', scientific(pvalue, digits = 1)))) +
  scale_fill_manual(values = colme) +
  labs(y = 'Bray-Curtis dissimilarity', x = '', fill = '', linetype = 'Phylum') +
  theme(legend.position = 'bottom', axis.ticks.x = element_blank()) 
bray_boxplot

b_time <- time_bray %>%
  ggplot() +
  geom_point(mapping = aes(x = as.numeric(date_dist, units = "days"), y = median, color = biota)) +
  geom_smooth(mapping = aes(x = as.numeric(date_dist, units = "days"), y = median, color = biota), 
              method = 'lm', se = TRUE, alpha = .2) +
  geom_text(bray_corr_time %>% 
              filter(biota == 'Bulk microbiota'), mapping = aes(x = 125, y = 0.68, 
                                                                label = paste('p=', round(p_value, digits = 2), '\n',
                                                                               'Corr:', round(corr, digits = 2)), color = biota)) +
  geom_text(bray_corr_time %>% 
              filter(biota == 'Ethanol treated sample'), mapping = aes(x = 165, y = 0.68, 
                                                                label = paste('p=', round(p_value, digits = 2), '\n',
                                                                              'Corr:', round(corr, digits = 2)), color = biota)) +
  
  scale_color_manual(values = colme) +
  labs(x = 'Days between sampling', y = 'Bray-Curtis dissimilarity', color = '') +
  theme(legend.position = 'bottom')
b_time

ggarrange(bray_boxplot + labs(tag='A'), b_time + labs(tag = 'B'), 
          common.legend = TRUE, legend = 'bottom', 
          widths = c(0.4, 0.6))

ggsave('longitudinal_amplicons/plots/EM_bray_box_time.png', dpi=600)

# Jaccard 
jaccard <- calc_dist(otutab, 'jaccard', 'M') %>%
  mutate(biota = 'Bulk microbiota') %>%
  rbind(calc_dist(otutab, 'jaccard', 'S') %>% 
          mutate(biota = 'Ethanol treated sample'))

# Kruskal test for Jaccard distances 
# Within individual 
jaccard_within <- filter(jaccard, same_person == 'Intra individual')
kruskal_within <- kruskal.test(value~biota, data = jaccard_within)
# But between which groups ?
wilcox_within <- pairwise.wilcox.test(jaccard_within$value, jaccard_within$biota, paired = FALSE)

# Between 
jaccard_between <- filter(jaccard, same_person == 'Inter individuals')
kruskal_between <- kruskal.test(value ~ biota, data = jaccard_between)
wilcox_between <- pairwise.wilcox.test(jaccard_between$value, jaccard_between$biota, paired = FALSE)

wilcox_jaccard <- data.frame(pvalue = wilcox_within$p.value, same_person = 'Intra individual') %>%
  rbind(data.frame(pvalue = wilcox_between$p.value, same_person = 'Inter individuals'))
colnames(wilcox_jaccard) <- c('pvalue', 'same_person')

time_jaccard <- jaccard %>%
  # Filter different individuals
  filter(same_person == 'Intra individual') %>%
  # group by difference between days and person
  group_by(biota, date_dist) %>%
  reframe(median=median(value)) %>%
  ungroup()

# Calculate correaltions between diff (time between samplings) and distance metric
jaccard_corr_time <- time_corr(time_jaccard)
jaccard_corr_time

jaccard_boxplot <- ggplot(jaccard) +
  geom_boxplot(mapping = aes(x = same_person, y = value, fill = biota)) +
  geom_text(wilcox_jaccard, mapping = aes( x= same_person, y = 1, label = paste('p=', scientific(pvalue, digits = 1)))) +
  scale_fill_manual(values = colme) +
  labs(y = 'Jaccard distance', x = '', fill = '') +
  theme(legend.position = 'bottom', axis.ticks.x = element_blank()) 
jaccard_boxplot

j_time <- time_jaccard %>%
  ggplot() +
  geom_point(mapping = aes(x = as.numeric(date_dist, units = "days"), y = median, color = biota)) +
  geom_smooth(mapping = aes(x = as.numeric(date_dist, units = "days"), y = median, color = biota), 
              method = 'lm', se = TRUE, alpha = .2) +
  geom_text(jaccard_corr_time %>% 
              filter(biota == 'Bulk microbiota'), mapping = aes(x = 125, y = 0.8, 
                                                                label = paste('p=', round(p_value, digits = 2), '\n',
                                                                              'Corr:', round(corr, digits = 2)), color = biota)) +
  geom_text(jaccard_corr_time %>% 
              filter(biota == 'Ethanol treated sample'), mapping = aes(x = 165, y = 0.8, 
                                                                       label = paste('p=', round(p_value, digits = 2), '\n',
                                                                                     'Corr:', round(corr, digits = 2)), color = biota)) +
  
  scale_color_manual(values = colme) +
  labs(x = 'Days between sampling', y = 'Jaccard distance', color = '') +
  theme(legend.position = 'bottom')
j_time

ggarrange(jaccard_boxplot + labs(tag='A'), j_time + labs(tag = 'B'), 
          common.legend = TRUE, legend = 'bottom', 
          widths = c(0.4, 0.6))

ggsave('longitudinal_amplicons/plots/EM_jaccard_box_time.png', dpi=600)

# NMDS plot
plot_nmds <- function(tab, method, i ) {
  # Separate the smaples 
  tab_filt <- as.data.frame(tab) %>%
    rownames_to_column('Group') %>%
    filter(substr(Group, 1, 1) == i) %>%
    column_to_rownames('Group') %>%
    as.matrix()
  
  # Calculate dissimilarity / distance
  dist <- vegdist(tab_filt, method = method)
  #  Calculate distances in 2D space
  nmds <- metaMDS(dist)
  # append metadata
  nmds_positions <-as.data.frame(scores(nmds, display='sites')) %>%
    rownames_to_column('Group') %>%
    left_join(metadata %>% select(Group, person, date), by = 'Group')
  return(nmds_positions)
}

nmds_bray <- plot_nmds(otutab, 'bray', 'M') %>%
  mutate(biota = 'Bulk microbiota') %>%
  rbind(plot_nmds(otutab, 'bray', 'S') %>%
          mutate(biota = 'Ethanol treated sample') )

ggarrange(plot_nmds(otutab, 'bray', 'M') %>%
  ggplot(aes(x= NMDS1, y= NMDS2, color = person)) +
  geom_point(size=4) +
  scale_x_continuous(breaks = c(-1, 0, 1)) +
  scale_y_continuous(breaks = c(-1, 0, 1)) +
  labs(x='', y='', color = 'Individual', subtitle = 'Bulk microbiota'), 
  plot_nmds(otutab, 'bray', 'S') %>%
    ggplot(aes(x= NMDS1, y= NMDS2, color = person)) +
    geom_point(size=4) +
    scale_x_continuous(breaks = c(-1, 0, 1)) +
    scale_y_continuous(breaks = c(-1, 0, 1)) +
    labs(x='', y='', color = 'Individual', subtitle = 'Ethanol treated samples'), 
  common.legend = TRUE, legend = 'bottom' )
ggsave('longitudinal_amplicons/plots/EM_nmds.png', dpi=600)

# If I normalize distances of each individual with min-max normalization, so that the dispersion of each individuals cluster does not account
# for the differences between microbiota and sporobiota!
# Normalized distances of each individual

# Bray-Curtis
dist_bray_norm = bray %>%
  filter(same_person == 'Intra individual') %>%
  group_by(person.x, biota) %>%
  # z-score normalization
  mutate(z_norm_value = ((value-mean(value)/sd(value))),
         # max-min normalization
         min_max_norm = (value - min(value))/(max(value) - min(value)))

z_score_plot = dist_bray_norm %>%
  ggplot(aes(x=biota, y=z_norm_value, fill=biota)) +
  geom_boxplot() +
  annotate("text", x=1.5, y= -1, label=paste("p=", scientific(wilcox.test(z_norm_value  ~ biota, data = dist_bray_norm)$p.value, digits =1)), size=3, color='black') +
  scale_fill_manual(values = colme) +
  labs(x='', y='Bray-Curtis dissimilarity (Z-score normalized)', fill='') +
  theme(axis.text.x = element_blank())

minmax_plot = dist_bray_norm %>%
  ggplot(aes(x=biota, y=min_max_norm, fill=biota)) +
  geom_boxplot() +
  scale_fill_manual(values = colme) +
  annotate("text", x=1.5, y= 1, label=paste("p=", scientific(wilcox.test(min_max_norm  ~ biota, data = dist_bray_norm)$p.value, digits =1)), size=3, color='black') +
  labs(x='', y='Bray-Curtis dissimilarity (min-max normalized)', fill='') +
  theme(axis.text.x = element_blank())

ggarrange(z_score_plot, minmax_plot,
          common.legend = TRUE,
          legend = 'bottom')
ggsave('longitudinal_amplicons/plots/EM_bray_norm_boxplot.png', dpi = 600)

# Jaccard 
dist_jaccard_norm = jaccard %>%
  filter(same_person == 'Intra individual') %>%
  group_by(person.x, biota) %>%
  # z-score normalization
  mutate(z_norm_value = ((value-mean(value)/sd(value))),
         # max-min normalization
         min_max_norm = (value - min(value))/(max(value) - min(value)))

z_score_plot = dist_jaccard_norm %>%
  ggplot(aes(x=biota, y=z_norm_value, fill=biota)) +
  geom_boxplot() +
  annotate("text", x=1.5, y= -1, label=paste("p=", scientific(wilcox.test(z_norm_value  ~ biota, data = dist_jaccard_norm)$p.value, digits =1)), size=3, color='black') +
  scale_fill_manual(values = colme) +
  labs(x='', y='Jaccard distances (Z-score normalized)', fill='') +
  theme(axis.text.x = element_blank())

minmax_plot = dist_jaccard_norm %>%
  ggplot(aes(x=biota, y=min_max_norm, fill=biota)) +
  geom_boxplot() +
  scale_fill_manual(values = colme) +
  annotate("text", x=1.5, y= 1, label=paste("p=", scientific(wilcox.test(min_max_norm  ~ biota, data = dist_jaccard_norm)$p.value, digits =1)), size=3, color='black') +
  labs(x='', y='Jaccard distances (min-max normalized)', fill='') +
  theme(axis.text.x = element_blank())

ggarrange(z_score_plot, minmax_plot,
          common.legend = TRUE,
          legend = 'bottom')
ggsave('longitudinal_amplicons/plots/EM_jaccard_norm_boxplot.png', dpi = 600)


# Additional test for usage of beta diveristy metrics: Is the difference we see between ethanol treated samples and bulk microbiota true, or 
# becouse of different relative abundances of OTUs (Bacillota have higher relative abundance)
sf <- otu_all %>%
  group_by(name) %>%
  reframe(mean_rel_abund =  mean(rel_abund), 
          sumsq_diff_abund = sum((outer(rel_abund, rel_abund, `-`)^2)[lower.tri(outer(rel_abund, rel_abund))])) %>%
  left_join(select(otu_all, name, biota), by = 'name')

ggplot(sf, aes(x = log10(mean_rel_abund), y = log10(sumsq_diff_abund), color = biota), alpha=.2) +
  geom_point() +
  scale_color_manual(values = colme) +
  labs(x = 'Mean relative abundance of OTU', y = 'Sum of squared differences between realtive abudnances of an OTU in different samples', color = '') +
  theme(legend.position = 'bottom')
ggsave('longitudinal_amplicons/plots/supplement_figure8.png', dpi = 600)

# Core community analysis 
# Core OTUs are those present in at least 11 time-points, irregardles of their relative or absolute abundance
core_otus = otu_all %>%
  mutate(PA = ifelse(value > 0, 1, 0)) %>%
  group_by(biota, person, name) %>%
  arrange(day, .by_group = TRUE) %>%
  # Create new column otu_sum is 1 if the OTU is present (PA > 0) on the current day and was not present on any of the previous days
  mutate(otu_cumsum = cumsum(PA)) %>%
  ungroup() %>%
  filter(otu_cumsum %in% c(11, 12)) %>%
  group_by(biota, person) %>%
  summarise(name = list(unique(name)), .groups = 'drop')

unnest(core_otus, name) %>%
  left_join(taxtab, by = 'name') %>%
  group_by(biota, person, Class) %>%
  summarise(number = n_distinct(name), .groups = 'drop') %>%
  mutate(class_fin = ifelse(number > 5, Class, 'Less than 5 OTUs per Class')) %>%
  ggplot(aes(x = biota, y = number, fill = class_fin)) +
  geom_bar(stat = 'identity') +
  facet_wrap(vars(person), scales = 'free_y') +
  labs( x = '', y= 'Number of unique OTUs', fill = 'Class')

# OTUs shared between all individuals 
core_all = unnest(core_otus, name) %>%
  mutate(value = 1) %>%
  pivot_wider(names_from = 'person', values_from = 'value', values_fill = 0) %>%
  group_by(biota, name) %>%
  summarise(sum_all = sum(A+B+C+D+E+F+G+H+I)) %>%
  filter(sum_all == 9) %>%
  left_join(taxtab, by = 'name') 

core_all %>%
  ggplot(aes(x = biota, fill = Class)) +
  geom_bar( stat = 'count') +
  labs(x = '', y = 'Number of OTUs present in all individuals across all time points')
ggsave('longitudinal_amplicons/plots/EMcoreotus_across_allPeople.png', dpi=600)


####
# Persistence of OTUs 
####

# How many OTUs are present in 1, 2, 3, 4 .. 12 samples and are unique? 
otu_all <- otu_all %>%
  mutate(PA = ifelse(value > 0, 1, 0))

fin = data.frame()
for (i in unique(otu_all$person)) {
  for (j in unique(otu_all$biota)) {
    # Select only one biota and person
    merged_sub = filter(otu_all, biota == j & person == i)
    merged_sub2 = select(merged_sub, name, Group, PA) %>%
      pivot_wider(names_from = 'name', values_from = 'PA') %>%
      column_to_rownames('Group') %>%
      t() %>%
      as.data.frame()
    # Calculate prevalence
    merged_sub2$prevalence = rowSums(merged_sub2)
    # Add metadata
    merged_sub3 = mutate(merged_sub2, 
                         biota = j, 
                         person = i, 
                         name = rownames(merged_sub2))
    rel_person <- otu_all %>%
      group_by(name, person, biota) %>%
      reframe(rel_abund = mean(rel_abund))
    
    merged_fin = left_join(select(merged_sub3, name, biota, person, prevalence), rel_person, by = c('name', 'person', 'biota'))
    fin = rbind(fin, merged_fin)
    
  }
}

fin_fin = fin %>% group_by(biota, person, prevalence) %>%
  filter(prevalence != 0) %>%
  summarise(count = sum(n_distinct(name)),
            names = list(unique(name)), 
            .groups = 'drop') 

fin_mean = fin_fin %>%
  group_by(biota, prevalence) %>%
  summarize(mean = mean(count),
            sd = sd(count), .groups = 'drop')

fin_percent = fin_fin %>% group_by(biota, person) %>%
  summarise(all = sum(count), .groups = 'drop') %>%
  full_join(fin_fin, by = join_by('biota', 'person')) %>%
  mutate(percent = (count/all) * 100)

fin_mean_percent = fin_percent %>%
  group_by(biota, prevalence) %>%
  summarize(mean = mean(percent),
            sd = sd(percent), .groups = 'drop')

ggplot() +
  geom_point(fin_percent, mapping = aes(x = prevalence, y = count, color=biota), size=3) +
  geom_line(fin_mean, mapping=aes(y=mean, x=prevalence, color = biota), linewidth=1.5) +
  scale_color_manual(values = colme) +
  scale_x_continuous(breaks = seq(0,14, by=1)) +
  labs(x='Occupancy (time-points of individual)', y= 'Number of OTUs present in # time-points', color = '') +
  theme(legend.position = 'bottom')
ggsave('longitudinal_amplicons/plots/EMoccupancy_count.png', dpi=600)

ggplot() +
  geom_point(fin_percent, mapping = aes(x = prevalence, y = percent, color = biota), size=3) +
  geom_line(fin_mean_percent, mapping=aes(y=mean, x=prevalence, color = biota), linewidth=1.5) +
  scale_color_manual(values = colme) +
  scale_x_continuous(breaks = seq(0,14, by=1)) +
  scale_y_continuous(breaks = seq(0,100, by=20)) +
  labs(x='Occupancy (time-points of individual)', y= 'Percent of OTUs present in # time-points (%)', color = '') +
  theme(legend.position = 'bottom')
ggsave('longitudinal_amplicons/plots/EM_occupancy_percent.png', dpi=600)

# How many OTUs are shared between people and how many are unique to only 1 person? 
core_otus <- otu_all_long %>%
  mutate(Group_clean = str_remove(Group, "-.*$")) %>%
  left_join(select(metadata, Group, date), by = join_by('Group_clean' == 'Group')) %>%
  group_by(name, person, fraction) %>%
  arrange(date, .by_group = TRUE) %>%
  # Create new column otu_sum is 1 if the OTU is present (PA > 0) on the current day and was not present on any of the previous days
  reframe(otu_cumsum = cumsum(PA)) %>%
  filter(otu_cumsum >= 4) %>%
  group_by(person, fraction) %>%
  summarise(name = list(unique(name)), .groups = 'drop') %>%
  mutate(core = sapply(name, function(x) length(unique(x)))) %>%
  left_join(all_person, by = 'person') %>%
  left_join(all_fraction, by = c('person', 'fraction'))

ggplot(core_otus, aes(x = person, y = (all_fraction/all_person) * 100, fill = fraction)) +
  geom_col() +
  labs(x = 'Individual', y = 'Fractions of microbiota (%)', fill = '')
ggsave('out/exploration/percent_fractions.png', width = 10, height = 20, unit = 'cm', dpi= 600)

ggplot(core_otus, aes(x = person, y = (core/all_person) * 100, fill = fraction)) +
  geom_col(position = 'dodge') +
  labs( x = '', y= 'OTUs in each fraction present in individual all the time (%)', fill = '')
ggsave('out/exploration/percent_core_fractions.png', width = 20, height = 15, units = 'cm', dpi=600)
