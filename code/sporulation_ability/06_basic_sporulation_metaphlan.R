# Sporulation abiliyty differentiate bacteria analysis 
# Library
library(dplyr)
library(readr)
library(tidyr)
library(ggplot2)
library(stringr)
library(tibble)
library(purrr)
library(ggpubr)
library(vegan)

set.seed(96)
theme_set(theme_bw(base_size = 14))

col <- c('#3CB371', '#f0a336')
col2 <- c('#A7E2C1', '#F7CD92')
colm <- '#3CB371'
cole <- '#f0a336'

# Sporulation ability 
# Before this part run analysis in the folder code/sporulation_ability
sporulation_ability <- read.table('data/sporulation_ability/sporulation_ability2021.tsv', sep = '\t', header = TRUE) %>% 
  as_tibble()

# Info on ethanol resistancy 
etoh_species <- read.table('data/longitudinal_shotgun/ethanol_resistant_SGB.tsv', sep = '\t', header = T)

metadata <- read.table('../longitudinal_shotgun/data/metadata.csv', sep = ';', header = TRUE) %>%  
  mutate(date = lubridate::dmy(date))


abund <- read_tsv('~/projects/longitudinal_shotgun/data/metaphlan_abundance_table.txt', comment = '#') %>%
  rename_with(~ str_remove(., '^profiled_'), starts_with('profiled_')) %>% 
  pivot_longer(-clade_name) %>% 
  filter(grepl('s__', clade_name), !grepl('t__', clade_name)) 


abund2 <- read_tsv('~/projects/longitudinal_shotgun/data/metaphlan_abundance_table.txt', comment = '#') %>%
  rename_with(~ str_remove(., '^profiled_'), starts_with('profiled_')) %>%
  filter(grepl('s__', clade_name), !grepl('t__', clade_name)) %>% 
  left_join(select(sporulation_ability, n_genes, PA, sporulation_ability, clade_name), by = 'clade_name') %>% 
  pivot_longer(-c(clade_name, PA, n_genes, sporulation_ability)) %>% 
  #mutate(clade_name = str_remove_all(clade_name, '[a-zA-Z]__')) %>%
  separate(clade_name, into=c('Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species'),
           sep="\\|") %>% 
  mutate(Phylum = ifelse(Phylum == 'p__Firmicutes', 'p__Bacillota', Phylum), 
         Domain = str_remove_all(Domain, 'k__'), 
         Phylum = str_remove_all(Phylum, 'p__'), 
         Class = str_remove_all(Class, 'c__'), 
         Order = str_remove_all(Order, 'o__'), 
         Family = str_remove_all(Family, 'f__'), 
         Genus = str_remove_all(Genus, 'g__'), 
         Species = str_remove_all(Species, 's__')) %>% 
  filter(name != 'MC013') %>% 
  left_join(metadata, by = join_by('name' == 'Group')) 


# What percentage of bacteria in the human gut can sporulate 
pre_abund <- abund2 %>% left_join(etoh_species, by = 'Species', relationship = 'many-to-many') %>% distinct()

pre_abund %>% 
  ggplot(aes(x = sporulation_ability, y = value)) +
  geom_boxplot() +
  scale_y_log10() +
  labs(x = '', y = 'Relative abundance [%]')
ggsave('out/sporulation/rel_abund_sporulation.png', dpi = 600)

pre_abund %>% 
  filter(biota == 'bulk microbiota', !is.na(sporulation_ability)) %>% 
  group_by(person, name, sporulation_ability, is_etoh_resistant) %>%  
  reframe(sum_rel = sum(value)) %>% 
  group_by(sporulation_ability, is_etoh_resistant) %>%
  reframe(mean_rel = mean(sum_rel), 
          sd_rel = sd(sum_rel))
  ggplot(aes(x = name, y = sum_rel, fill = sporulation_ability)) +
  geom_col() +
  facet_wrap(is_etoh_resistant~person, scales = 'free_x') +
  theme(axis.text.x = element_text(angle = 90))

  # Number of sporeformers in spore form 
pre_abund %>% 
  filter(biota == 'bulk microbiota', !is.na(sporulation_ability)) %>% 
  filter(value > 0) %>% 
  group_by(sporulation_ability, is_etoh_resistant) %>% 
  reframe(PA = n_distinct(Species)) 

# number of spperson# number of species 
pre_abund %>% 
  mutate(PA = ifelse(value > 0, 1, 0)) %>% 
  filter(PA == 1) %>% 
  group_by(sporulation_ability) %>% 
  reframe(PA = n_distinct(Species)) 

# Non-spore-former      292
# Spore-former          235
# NA                    506

abund3 <- filter(abund2, !is.na(sporulation_ability), Domain == 'Bacteria', biota == 'bulk microbiota') %>% 
  left_join(etoh_species, by = 'Species', relationship = 'many-to-many') %>% 
  mutate(is_etoh_resistant = dplyr::if_else(is_etoh_resistant == TRUE, 'Ethanol-resistant', 'Non-ethanol resistant', missing = 'Non-ethanol resistant')) %>%  
  distinct()

unique(abund3$is_etoh_resistant)

# relative abudnance of spore/ non-spore forming 
rel <- abund3 %>% 
  ggplot(aes(x = value, y = Phylum, fill = sporulation_ability)) +
  geom_boxplot() +
  scale_x_log10() +
  scale_fill_manual(values = col) +
  labs(x = 'Relative abundance [log10]', y = '', fill = '') +
  theme(legend.position = 'bottom')
rel
ggsave('out/sporulation/SNS_rel_abund.png')

abund4 <- abund3 %>%
  mutate(spore_etoh = case_when(
    sporulation_ability == 'Spore-former' & is_etoh_resistant ==  'Ethanol-resistant'  ~ 'Ethanol-resistant spore-former',
    sporulation_ability == 'Spore-former' & is_etoh_resistant == 'Non-ethanol resistant'  ~ 'Non-ethanol resistant spore-former',
    sporulation_ability == 'Non-spore-former' & is_etoh_resistant ==  'Ethanol-resistant' ~ 'Ethanol resistant non-sporeforming bacteria',
    TRUE ~ 'Non-ethanol resistant non-spore forming bacteria'))

unique(abund4$spore_etoh)

abund4 %>%
  ggplot(aes(x = value, y = Phylum, fill = spore_etoh)) +
  geom_boxplot() +
  scale_x_log10() +
  #scale_fill_manual(values = col) +
  labs(x = 'Relative abundance [log10]', y = '', fill = '') +
  theme(legend.position = 'bottom') +
  guides(fill = guide_legend(nrow = 2))
  
abund3 %>% 
  ggplot(aes(x = value, y = Phylum, fill = sporulation_ability)) +
  geom_boxplot() +
  scale_x_log10() +
  scale_fill_manual(values = col) +
  labs(x = 'Relative abundance [log10]', y = '', fill = '') +
  facet_wrap(~is_etoh_resistant, scales = 'free') +
  theme(legend.position = 'bottom') 
ggsave('out/sporulation/rel_etoh_sporulation_Phylum.png', dpi = 600)

abund3 %>% 
  ggplot(aes(x = value, y = Phylum, fill = is_etoh_resistant)) +
  geom_boxplot() +
  scale_x_log10() +
  #scale_fill_manual(values = col) +
  labs(x = 'Relative abundance [log10]', y = '', fill = '') +
  theme(legend.position = 'bottom') +
  stat_compare_means()
ggsave('out/sporulation/etoh_relabund.png')  

# Is distribution of relative abundances for Bacillota the same for non-spore and spore-forming? 
bacillota <- filter(abund3, Phylum == 'Bacillota')

wilcox.test(filter(bacillota, sporulation_ability == 'Spore-former')$value, filter(bacillota, sporulation_ability == 'Non-spore-former')$value, alternative = 'greater')

# Wilcoxon rank sum test with continuity correction
# 
# data:  filter(rel_bacillota, sporulation_ability == "Spore-former")$value and filter(rel_bacillota, sporulation_ability == "Non-spore-former")$value
# W = 257906560, p-value < 2.2e-16
# alternative hypothesis: true location shift is greater than 0

# Number of species 
no <- abund3 %>% 
  mutate(PA = ifelse(value > 0, 1, 0)) %>%  
  filter(PA == 1) %>% 
  group_by(sporulation_ability, Phylum) %>% 
  reframe(sum = n_distinct(Species)) %>% 
  ggplot(aes(x = sum, y = Phylum, fill = sporulation_ability)) +
  geom_col(position = position_dodge(width = 0.9)) +
  geom_text(aes(label = sum),
            position = position_dodge(width = 0.9), hjust = -0.1) +
  scale_fill_manual(values = col) +
  labs(x = '# Species', y = '', fill = '') +
  theme_bw(base_size=14) +
  theme(legend.position = 'bottom') 
  
no
ggsave('out/sporulation/SNS_N_species.png')

# density rel abund 
bacillota %>% 
  ggplot(aes(x = value, color = sporulation_ability)) +
  geom_density(linewidth = 2) +
  scale_x_log10() +
  labs(x = 'Relative abundance [log10]', y = 'Density', color = '')
ggsave('out/sporulation/SNS_density_relabund_bacillota.png')


ggarrange(no + labs(tag = 'A'), 
          rel + labs(tag = 'B') + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()), 
          common.legend = T, legend = 'bottom', widths = c(1, .6))
ggsave('out/sporulation/number_species_rel_abund.png', dpi = 600)

# Prevalence 
preval3 <- abund3 %>% 
  group_by(sporulation_ability, person, Phylum, Species) %>% 
  reframe(all_timepoints = n(), 
          timepoints_present = sum(value > 0),
          prevalence = (timepoints_present / all_timepoints) * 100) %>%
  # Calculate number of OTUs per person x treatment x Phylum
  filter(prevalence > 0) %>% 
  group_by(person) %>%
  mutate(all_species_per_person = n_distinct(Species)) %>%
  ungroup() %>%
  group_by(person, sporulation_ability, Phylum, prevalence) %>%
  reframe(n_species_person_spore_phylum = n_distinct(Species), 
          per_species = (n_species_person_spore_phylum / all_species_per_person) * 100) %>% 
  distinct() 

# Statistics for prevalence ? 
wt <- abund3 %>% 
  group_by(sporulation_ability, person, Species) %>% 
  reframe(all_timepoints = n(), 
          timepoints_present = sum(value > 0),
          timepoints_missing = sum(value = 0),
          prevalence = (timepoints_present / all_timepoints) * 100) 

wt_results <- wilcox.test(filter(wt, sporulation_ability == 'Spore-former')$prevalence, filter(wt, sporulation_ability == 'Non-spore-former')$prevalence)
wt_results
# Wilcoxon rank sum test with continuity correction
# 
# data:  filter(wt, sporulation_ability == "Spore-former")$prevalence and filter(wt, sporulation_ability == "Non-spore-former")$prevalence
# W = 3265983, p-value < 2.2e-16
# alternative hypothesis: true location shift is not equal to 0
# By person 
group_by(wt, person) %>% 
  reframe(wilcox.test(filter(wt, sporulation_ability == 'Spore-former')$prevalence, filter(wt, sporulation_ability == 'Non-spore-former')$prevalence)$p.value)

# 1 A          1.07e-31
# 2 B          1.07e-31
# 3 C          1.07e-31
# 4 D          1.07e-31
# 5 E          1.07e-31
# 6 F          1.07e-31
# 7 G          1.07e-31
# 8 H          1.07e-31
# 9 I          1.07e-31

preval_plot <- preval3 %>% 
  filter(Phylum == 'Bacillota') %>% 
  ggplot(aes(x = prevalence, y = per_species, color = person)) +
  geom_line(linewidth = 2) +
  geom_point(size = 2) +
  #annotate("text", x = 50, y = 10, label = paste0('Wilcox test; p < ', scales::scientific(wt_results$p.value, digits = 3)), color = 'black') +
  facet_wrap(~sporulation_ability, nrow = 1, scales = 'free_y') +
  #scale_color_manual(values = col2) +
  labs(x = 'Prevalence [% days present within person]', y = expression(paste(italic("Bacillota"), " species [% of species at this prevalence]")), color = '') +
  theme(legend.position = 'bottom') +
  guides(color = guide_legend(nrow = 1))
preval_plot
ggsave('out/sporulation/SNS_prevalence_person_Bacillota.png')

preval3 %>% 
  filter(Phylum == 'Bacillota') %>%
  ggplot(aes(x = prevalence, y = n_species_person_spore_phylum, color = person)) +
  geom_point(size = 3) +
  #annotate("text", x = 50, y = 10, label = paste0('Wilcox test; p < ', scales::scientific(wt_results$p.value, digits = 3)), color = 'black') +
  facet_wrap(~sporulation_ability, nrow = 1, scales = 'free_y') +
  #scale_color_manual(values = col2) +
  labs(x = 'Prevalence [% days present within person]', y = expression(paste('#', italic("Bacillota"), " species")), color = '') +
  theme(legend.position = 'bottom') +
  guides(color = guide_legend(nrow = 1))
ggsave('out/sporulation/SNS_prevalence_number_Bacillota.png')

preval3 %>% 
  ggplot(aes(x = prevalence, y = n_species_person_spore_phylum, color = person)) +
  geom_point(size = 3) +
  #annotate("text", x = 50, y = 10, label = paste0('Wilcox test; p < ', scales::scientific(wt_results$p.value, digits = 3)), color = 'black') +
  facet_wrap(~sporulation_ability, nrow = 1, scales = 'free_y') +
  #scale_color_manual(values = col2) +
  labs(x = 'Prevalence [% days present within person]', y = "# species", color = '') +
  theme(legend.position = 'bottom') +
  guides(color = guide_legend(nrow = 1))
ggsave('out/sporulation/SNS_n_species_prevalence_person.png')

###
# Alpha diveristy 
#bacillota <- filter(abund3, Phylum == 'Bacillota')
bacillota <- mutate(abund3, name = paste0(ifelse(sporulation_ability == 'Spore-former', 'S_', 'NS_'),name))

n <- filter(bacillota, value > 0) %>% 
  group_by(name) %>% 
  reframe(richness = n_distinct(Species)) 

tab <- bacillota %>% 
  group_by(Species, name) %>%  
  reframe(value = mean(value)) %>% 
  pivot_wider(names_from = 'name', values_from = 'value', values_fill = 0) %>% 
  column_to_rownames('Species') %>% 
  t()

shannon = diversity(tab, index = 'shannon')

alpha <- left_join(n, as_tibble(as.list(shannon)) %>% 
                     pivot_longer(names_to = 'name', values_to = 'shannon', cols = starts_with(c('NS_', 'S_'))), 
                   by = 'name') %>%
  left_join(bacillota, by = 'name') %>% 
  mutate(person2 = person) 

event_data <- read.table('data/extreme_event_data.csv', sep = ',', header = TRUE)

# richness
ggplot(alpha, aes(x = day, y = richness)) +
  geom_rect(data = event_data, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = extremevent_type), inherit.aes = FALSE,
            alpha = 0.6) +
  scale_fill_manual(values = c('white','#d94343', '#d98e43', '#f1f011', '#0c9910', '#3472b7', '#7934b7', '#b73485', '#0f5618')) +
  geom_line(data = dplyr::select(alpha, -person) %>% filter(sporulation_ability == 'Non-spore-former'), 
            aes(group=person2), color= colm, linewidth=0.5, alpha=0.5) +
  geom_line(data=alpha %>% dplyr::select(-person) %>% filter(sporulation_ability == 'Spore-former'), 
            aes(group=person2), color= cole, linewidth=0.5, alpha=0.5) +
  geom_line(data=alpha %>% filter(sporulation_ability == 'Non-spore-former'),
            aes(color=person), color= colm, linewidth=1.2) +
  geom_line(data=alpha %>% filter(sporulation_ability == 'Spore-former'), 
            color=cole, linewidth=1.2) +
  facet_wrap(~person, scales = 'free') +
  labs(x='Day', y= 'Richness', fill = 'Event')
ggsave('out/sporulation/SNS_richness_all.png')

ggplot(alpha, aes(x=day, y=shannon)) +
  geom_rect(data = event_data, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = extremevent_type), inherit.aes = FALSE,
            alpha = 0.6) +
  scale_fill_manual(values = c('white','#d94343', '#d98e43', '#f1f011', '#0c9910', '#3472b7', '#7934b7', '#b73485', '#0f5618')) +
  geom_line(data = dplyr::select(alpha, -person) %>% filter(sporulation_ability == 'Non-spore-former'), 
            aes(group=person2), color= colm, linewidth=0.5, alpha=0.5) +
  geom_line(data=alpha %>% dplyr::select(-person) %>% filter(sporulation_ability == 'Spore-former'), 
            aes(group=person2), color= cole, linewidth=0.5, alpha=0.5) +
  geom_line(data=alpha %>% filter(sporulation_ability == 'Non-spore-former'),
            aes(color=person), color= colm, linewidth=1.2) +
  geom_line(data=alpha %>% filter(sporulation_ability == 'Spore-former'), 
            color=cole, linewidth=1.2) +
  facet_wrap(~person, scales = 'free') +
  labs(x='Day', y= 'Shannon', fill = 'Event')
ggsave('out/sporulation/SNS_shannon_all.png')

# Beta diversity & density of clustering
dist <- vegdist(tab, method = 'bray')
mds <- metaMDS(dist)

mds_data <- as.data.frame(mds$points) %>%  
  rownames_to_column('name') %>% 
  mutate(name_old = str_remove_all(name, 'S_'), 
         name_old = str_remove_all(name_old, 'N')) %>% 
  left_join(metadata, by = join_by('name_old' == 'Group')) %>% 
  mutate(spore_ability = substr(name, 1, 2))

ggplot(mds_data, aes(x = MDS1, y = MDS2, color = person, shape = spore_ability)) +
  geom_point(size = 3)

# 
dist_bray <- as.matrix(dist) %>%
  as_tibble(rownames = 'Group') %>%
  pivot_longer(-Group) %>%
  filter(Group != name) %>%  
  mutate(sample_pairs = paste(Group, name)) %>%
  group_by(sample_pairs) %>%
  reframe(mean_value = mean(value, na.rm = TRUE), 
          median_value = median(value, na.rm = TRUE)) 

# Tidy the Bray data and join with metadata
bray <- dist_bray %>%
  separate(sample_pairs, into = c("Group", "name"), sep = " ") %>%
  mutate(Group_clean = str_remove(Group, "^(NS_|S_)"),
         name_clean = str_remove(name, "^(NS_|S_)"), 
         spore_ability.x = substr(Group, 1, 2), 
         spore_ability.y = substr(name, 1, 2)) %>%
  left_join(metadata %>% select(Group, person, date), by = join_by('Group_clean' == 'Group')) %>%
  left_join(metadata %>% select(Group, person, date), by = join_by('name_clean' == 'Group')) %>%
  mutate(same_person = ifelse(person.x == person.y, 'Within individual', 'Between individuals'),
         same_fraction = ifelse(spore_ability.x == spore_ability.y, 'Yes', 'No')) %>%
  filter(same_fraction == 'Yes')

one <- bray %>%
  mutate(spore_ability.x = ifelse(spore_ability.x == 'S_', 'Spore-forming', 'Non-spore-forming')) %>% 
  ggplot(aes(x=spore_ability.x, y=median_value, fill=spore_ability.x)) +
  geom_boxplot() +
  stat_compare_means() +
  scale_fill_manual(values = col) +
  facet_wrap(~same_person) +
  labs(x = '', y = 'Median Bray-Curtis dissimilarity', fill = '') 
one 
ggsave('out/sporulation/SNS_bray_curtis_boxplot_all.png')

# in time 
time_bray <-  as.matrix(dist) %>%
  as_tibble(rownames = 'Group') %>%
  pivot_longer(-Group) %>%
  filter(Group != name) %>%  
  mutate(Group_clean = str_remove(Group, "^(NS_|S_)"),
         name_clean = str_remove(name, "^(NS_|S_)"), 
         spore_ability.x = substr(Group, 1, 2), 
         spore_ability.y = substr(name, 1, 2)) %>% 
  left_join(metadata %>% select(Group, person, date), by = join_by('Group_clean' == 'Group')) %>%
  left_join(metadata %>% select(Group, person, date), by = join_by('name_clean' == 'Group')) %>%
  mutate(same_person = ifelse(person.x == person.y, 'Within individual', 'Between individuals'),
         same_fraction = ifelse(spore_ability.x == spore_ability.y, 'Yes', 'No'),
         spore_ability.x = ifelse(spore_ability.x == 'S_', 'Spore-forming', 'Non-spore-forming')) %>%
  filter(same_fraction == 'Yes') %>%  
  # Filter different individuals
  filter(same_person == 'Within individual') %>%
  # Calculate the difference between sampling times
  mutate(diff=abs(date.x-date.y)) %>%
  # group by difference between days and person
  group_by(spore_ability.x, person.x, diff) %>%
  reframe(median=median(value)) 

time <- time_bray %>%
  ggplot(aes(x=diff, y=median, color=spore_ability.x)) +
  geom_point() +
  geom_smooth(method = 'lm') +
  stat_cor(method = 'pearson', aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~")), size = 4) +
  labs(x='Days between sampling', y='Median Bray-Curtis distance', color='') +
  theme(legend.position = 'bottom') +
  guides(fill = guide_legend(ncol = 2)) +
  scale_color_manual(values = col) 
time 
ggsave('out/sporulation/SNS_time_braycurtis_all.png')

ggarrange(one + labs(tag = 'A'), time + labs(tag = 'B'), common.legend = T,  
          legend = 'bottom')
ggsave('out/sporulation/intra_inter_time_BC.png', dpi = 600)

# Sharing species between individuals 
# not by Phylum 
n_species <- bacillota %>% 
  mutate(PA = ifelse(value > 0, 1, 0)) %>% 
  group_by(sporulation_ability, person, Species) %>%
  reframe(present = sum(PA == 1)) %>%  
  filter(present > 12) %>% 
  group_by(sporulation_ability) %>% 
  reframe(n = n_distinct(Species))

n_species_people <- bacillota %>% 
  mutate(PA = ifelse(value > 0, 1, 0)) %>% 
  group_by(sporulation_ability, person, Species) %>%
  reframe(present = sum(PA == 1)) %>%  
  filter(present > 12) %>% 
  mutate(present = ifelse(present > 0, 1, 0)) %>% 
  select(Species, person, present, sporulation_ability) %>% 
  pivot_wider(names_from = 'person', values_from = 'present', values_fill = 0) %>%  
  mutate(n_people = A+B+C+D+E+F+G+H+I, 
         n_people = ifelse(n_people == 1, 'Present in 1 individual', 'Present > 1 individual')) %>% 
  group_by(sporulation_ability, n_people) %>% 
  reframe(n_species = n_distinct(Species)) %>% 
  left_join(n_species, by = 'sporulation_ability') %>% 
  mutate(per_species = (n_species/n)*100) 

n_species_people %>% 
  #mutate(n_people = factor(n_people, levels = c('Present in 1 individual', 'Present > 1 individual'))) %>% 
  ggplot(aes(y = sporulation_ability, x = per_species, fill = n_people))+
  geom_col() +
  scale_fill_manual(values = c('#F2933F', '#3F9EF2')) +
  labs(y = '', x = 'OTUs [%]', fill = '') +
  theme(legend.position = 'bottom') +
  guides(fill=guide_legend(nrow=1,byrow=TRUE))
ggsave('out/sporulation/present_one_more_people.png', dpi = 600)


##
##
# Ethanol resistant and spore-forming = active spore-formers who are they? 


# How many spore-formers are in the ethanol-treated and untreated samples? 

abund5 <- read_tsv('~/projects/longitudinal_shotgun/data/metaphlan_abundance_table.txt', comment = '#') %>%
  rename_with(~ str_remove(., '^profiled_'), starts_with('profiled_')) %>%
  mutate(clade_name2 = clade_name) %>% 
  filter(grepl('s__', clade_name), !grepl('t__', clade_name)) %>% 
  left_join(select(sporulation_ability, n_genes, PA, sporulation_ability, clade_name), by = 'clade_name') %>% 
  pivot_longer(-c(clade_name, clade_name2, PA, n_genes, sporulation_ability)) %>% 
  #mutate(clade_name = str_remove_all(clade_name, '[a-zA-Z]__')) %>%
  separate(clade_name, into=c('Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species'),
           sep="\\|") %>% 
  mutate(Phylum = ifelse(Phylum == 'p__Firmicutes', 'p__Bacillota', Phylum), 
         Domain = str_remove_all(Domain, 'k__'), 
         Phylum = str_remove_all(Phylum, 'p__'), 
         Class = str_remove_all(Class, 'c__'), 
         Order = str_remove_all(Order, 'o__'), 
         Family = str_remove_all(Family, 'f__'), 
         Genus = str_remove_all(Genus, 'g__'), 
         Species = str_remove_all(Species, 's__')) %>% 
  filter(name != 'MC013') %>% 
  left_join(metadata, by = join_by('name' == 'Group')) 


abund5 %>% filter(sporulation_ability == 'Spore-former', value > 0) %>% 
  group_by(biota) %>%  
  reframe(n_spore = n_distinct(clade_name2))

abund5 %>% filter(sporulation_ability == 'Non-spore-former', value > 0) %>% 
  group_by(biota) %>%  
  reframe(n_spore = n_distinct(clade_name2))

abund5 %>%  filter(value > 0) %>% 
  group_by(biota) %>%  
  reframe(n = n_distinct(clade_name2))

# What about relative abundance 
abund5 %>% 
  group_by(name, biota, person, day, sporulation_ability) %>%  
  reframe(rel = sum(value)) %>%  
  group_by(biota, sporulation_ability) %>% 
  reframe(mean= mean(rel), 
          sd = sd(rel))

# Is there a statisticaly significant enrichment of spore-formers in the ethanol-treated samples?
n_rel <- abund5 %>% 
  mutate(sample = ifelse(biota == 'bulk microbiota', 'untreated sample', 'ethanol-treated sample')) %>%
  mutate(pa = ifelse(value > 0, 1, 0)) %>% 
  group_by(name, sporulation_ability, sample) %>% 
  reframe(sum_n = sum(pa), 
          sum_rel = sum(value)) %>% 
  filter(sporulation_ability == 'Spore-former') %>% 
  mutate(original_sample = substr(name, 2, 5)) %>% 
  select(-name)

n_rel2 <- filter(n_rel, sample == 'untreated sample') %>% 
  left_join(filter(n_rel, sample == 'ethanol-treated sample'), by = 'original_sample') %>% 
  filter(!is.na(sample.y))


# Number of species detected
wilcox.test(n_rel2$sum_n.x, n_rel2$sum_n.y, paired = TRUE, alternative = 'greater')
# Wilcoxon signed rank test with continuity correction
# 
# data:  n_rel2$sum_n.x and n_rel2$sum_n.y
# V = 3655, p-value = 5.913e-16
# alternative hypothesis: true location shift is greater than 0


# Relative abundance 
wilcox.test(n_rel2$sum_rel.x, n_rel2$sum_rel.y, paired = TRUE, alternative = 'greater') 
# Wilcoxon signed rank test with continuity correction
# 
# data:  n_rel2$sum_rel.x and n_rel2$sum_rel.y
# V = 1125, p-value = 0.999
# alternative hypothesis: true location shift is greater than 0
