library(dplyr)
library(readr)
library(tibble)
library(tidyr)
library(ggplot2)

set.seed(96)
theme_set(theme_bw())

# Metadata 
metadata = read.table('data/metadata.csv', header= TRUE, sep = ';') %>%
  mutate(date = dmy(date)) 

contigs_virueses <- readRDS('data/intermediate/contigs_viruses.RDS')
contigs <- readRDS('data/intermediate/contigs.RDS')

# Import file of Virus indentification and clasification results for viruses from program geNomad. 
# Viruses can be taxonomically assigned up to the family level, but not to specific genera or species within that family. 
# The taxonomy is presented with a fixed number of fields (corresponding to taxonomic ranks) separated by semicolons, with empty fields left blank.
g_viruses <- read_tsv('data/geNomad/01.sqmMicrobiota_virus_summary.tsv') %>%
  separate(taxonomy, into=c('Kingdom', 'Clade', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species'), sep=";")

# How likely is for this sequence to be really virus ?
g_viruses %>%
  ggplot(aes(x = virus_score)) +
  geom_bar(stat = 'count')
ggsave('plots/geNomad/virus_score.png')

g_viruses %>%
  ggplot(aes(y = reorder(seq_name, marker_enrichment), x = marker_enrichment)) +
  geom_point() +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
  labs(y = 'Contigs', x = 'total enrichment of viral markers in a contig')
ggsave('plots/geNomad/marker_enrichment.png')

# Filter contigs that had virus score more than 0.8 and had at least 1% of genes recognised as marker genes
virus <- filter(g_viruses, virus_score > 0.8 & marker_enrichment > 1) %>%
  left_join(contigs, by = join_by('seq_name' == 'contigID')) 


# How many contigs are determined as ... 
# Crucial plots for results od disertation
contigs_long <- contigs %>% 
  mutate(kingdom = ifelse(contigID %in% virus$seq_name, 'Viruses', Kingdom)) %>%
  select(Kingdom, kingdom, starts_with('Raw read count M')) %>%
  mutate(reads_across_samples = sum(c_across(starts_with('Raw read ')), na.rm = TRUE)) 

contigs_long_SQM <- contigs_long %>%
  group_by(Kingdom) %>%
  summarise(no_reads = sum(reads_across_samples), .groups = 'drop') %>%
  mutate(x='x', 
         per_reads = no_reads/sum(no_reads)*100) 

contigs_long_geNomad <- contigs_long %>%
  group_by(kingdom) %>%
  summarise(no_reads = sum(reads_across_samples), .groups = 'drop') %>%
  mutate(x='x', 
         per_reads = no_reads/sum(no_reads)*100) 


ggarrange(contigs_long_SQM %>% ggplot(aes(x = x, y = per_reads, fill = Kingdom, label = paste(round(per_reads, digits = 2), '%'))) +
            geom_col() +
            geom_text(position = position_stack(vjust = 0.5), size = 3.5) +
            scale_fill_manual(values = c('red', '#3da217', '#426dc9', '#e353d4', 'grey')) +
            labs(x = '', y = '% reads', fill = '', caption = 'DIAMOND with LCA') +
            theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()), 
          contigs_long_geNomad %>% ggplot(aes(x = x, y = per_reads, fill = kingdom, label = paste(round(per_reads, digits = 2), '%'))) +
            geom_col() +
            geom_text(position = position_stack(vjust = 0.5), size = 3.5) +
            scale_fill_manual(values = c('red', '#3da217', '#426dc9', '#e353d4', 'grey')) +
            labs(x = '', y = '% reads', fill = '', caption = 'DIAMOND with LCA & geNomad') +
            theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()), 
          abund %>% filter(is.na(Phylum)) %>%
            group_by(Kingdom) %>%
            mutate(mean_rel = mean(c_across(starts_with('M')), na.rm = TRUE), 
                   Kingdom = ifelse(Kingdom == 'UNCLASSIFIED', NA, Kingdom)) %>%
            select(Kingdom, Phylum,  mean_rel) %>%
            rbind(data.frame(Kingdom = c('Viruses'), mean_rel = 0, Phylum = NA)) %>%
            ggplot(aes(x = Phylum, y = mean_rel, fill = Kingdom, label = ifelse(round(mean_rel, digits = 2) > 0, paste(round(mean_rel, digits = 2), '%'), ''))) +
            geom_col() +
            geom_text(position = position_stack(vjust = 0.5), size = 3.5) +
            scale_fill_manual(values = c( 'red', '#3da217', '#426dc9', '#e353d4', 'grey')) +
            labs(x = '', y = 'Relative abundance', fill = '', caption = 'metaphlan') +
            theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()), 
          common.legend = TRUE, legend = 'bottom', nrow = 1)
ggsave('plots/contigs/rel_abund_kingdom_3.png')


# 
kingdom_match <- mutate(virus, kingdom_match = ifelse(Kingdom.x == Kingdom.y, 'match', '')) %>%
  filter(kingdom_match == 'match')

virus_reads_long <- virus %>%
  select(seq_name, Kingdom.x, Clade.x, Phylum.x, Class.x, Order.x, Family.x, 
         Kingdom.y, Clade.y, Phylum.y, Class.y, Order.y, Family.y, starts_with('Raw read')) %>%
  pivot_longer(names_to = 'Group', values_to = 'value', cols = starts_with('Raw read')) %>%
  mutate(Group = str_remove_all(Group, 'Raw read count '))

virus_reads_long %>%
  group_by(seq_name, Kingdom.x, Clade.x, Phylum.x, Class.x, Order.x, Family.x, Kingdom.y, Clade.y, Phylum.y, Class.y, Order.y, Family.y) %>%
  summarise(value = sum(value, na.rm = TRUE)) %>%
  pivot_longer(names_to = 'Kingdom', values_to = 'which', cols= c(Kingdom.x, Kingdom.y)) %>%
  mutate(which = ifelse(is.na(which), 'Unclassified', which), 
         Kingdom = ifelse(Kingdom == 'Kingdom.x', 'geNomad', 'SQM')) %>%
  ggplot(aes(x = Kingdom, y = value, fill = which)) +
  geom_col() +
  labs(fill = '', x = '', y = '# reads')
ggsave('plots/geNomad/geNomad_SQM.png')

# What kind of viruses do we have, by person!
virus_reads_long_meta <- virus_reads_long %>%
  left_join(select(metadata, Group, person, biota, date, time_point), by = 'Group') 

virus_reads_long_meta%>%
  ggplot(aes( x = time_point, y = value, fill = Order.x)) +
  geom_col() +
  facet_grid(biota ~ person)
ggsave('plots/geNomad/Order_person.png')


virus_reads_long_meta %>%
  filter(biota == 'bulk microbiota' & !is.na(Order.x) & Order.x != '' ) %>%
  group_by(person, time_point, date, Order.x) %>%
  summarise(sum_value = sum(value, na.rm = TRUE), .groups = 'drop') %>%
  ggplot(aes(x = time_point, y = sum_value, color = Order.x)) +
  geom_point(size=2) +
  geom_line(linewidth=1) +
  scale_y_log10() +
  facet_wrap(~ person, scales = 'free_y') +
  labs(x = 'Time point', y = 'log10(# reads)', color = 'Order')
ggsave('plots/geNomad/order_by_time_log.png')

alpha_viruses <- virus_reads_long_meta %>%
  filter(!is.na(Order.x) & Order.x != '' ) %>%
  group_by(Group, Order.x) %>%
  summarise(value = sum(value, na.rm = TRUE), .groups = 'drop') %>%
  pivot_wider(names_from = 'Order.x', values_from = 'value') %>%
  column_to_rownames('Group') %>%
  estimateR() %>%
  t() %>%
  as.data.frame() %>%
  rownames_to_column('Group') %>%
  left_join(alpha_meta, by = 'Group')

# Correlation of alpha diveristy of viruses / OTUs 

alpha_viruses %>%
  pivot_longer(names_to = 'which', values_to = 'Observed', cols = c('S.obs.x', 'S.obs.y')) %>%
  mutate(which = ifelse( which == 'S.obs.x', 'Order viruses', 'OTU')) %>%
  filter(!is.na(biota) & !is.na(person)) %>%
  ggplot(aes(x = time_point, y = Observed, colour = which)) +
  geom_point(size=2) +
  geom_line(linewidth=1) +
  scale_y_log10() +
  facet_grid(biota~ person) +
  labs(x = 'Time point', y = '# OTUs / order viruses', color = '')
ggsave('plots/geNomad/alpha_corr_log.png')

cor.test(filter(alpha_viruses, biota == 'bulk microbiota')$S.obs.x, filter(alpha_viruses, biota == 'bulk microbiota')$S.obs.y, method = 'pearson')

cor.test(filter(alpha_viruses, biota == 'ethanol treated sample')$S.obs.x, filter(alpha_viruses, biota == 'ethanol treated sample')$S.obs.y, method = 'pearson')

# ethanol tretaed look like it could be correlational, but what if I separate by individual 
corr_res <- data.frame()

for (i in unique(alpha_viruses$person)) {
  alpha_filt <- filter(alpha_viruses, person == i)
  correlation <- cor.test(filter(alpha_filt, biota == 'ethanol treated sample')$S.obs.x, 
                          filter(alpha_filt, biota == 'ethanol treated sample')$S.obs.y, 
                          method = 'pearson')
  corr_res <- rbind(corr_res, data.frame(person = i, 
                                         pvalue = correlation$p.value, 
                                         correlation = correlation$estimate))
}


corr_res


# Beta diversity
viruses_tab <- virus_reads_long %>%
  filter(!is.na(Order.x) & Order.x != '' ) %>%
  select(Order.x, value, Group) %>%
  group_by(Order.x, Group) %>%
  summarise(value = sum(value, na.rm = TRUE)) %>%
  pivot_wider(names_from = 'Order.x', values_from = 'value') %>%
  column_to_rownames('Group')

viruses_beta <- vegdist(viruses_tab, method = 'bray')

nmds <-  metaMDS(viruses_beta)
nmds_positions <- as.data.frame(scores(nmds, display='sites')) %>%
  rownames_to_column('Group') %>%
  left_join(metadata %>% select(Group, person, date, biota), by = 'Group')

nmds_positions %>%
  ggplot(aes(x=NMDS1, y=NMDS2, color=person, shape = biota)) +
  geom_jitter(size=4) +
  facet_grid(~biota, scales = 'free') +
  labs(x='', y='', color='Individual')
ggsave('plots/geNomad/beta_viruses.png')

