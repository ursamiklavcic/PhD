# geNomad
library(stringr)
library(ggplot2)
library(readr)
library(dplyr)
library(tibble)
library(lubridate)
library(tidyr)

theme_set(theme_bw(base_size=14))

metadata <- read.table('~/projects/longitudinal_shotgun/data/metadata.csv', header= TRUE, sep = ';') %>%
  mutate(date = dmy(date))

viruses <- read.table('~/projects/longitudinal_shotgun/data/geNomad/01.sqmMicrobiota_virus_summary.tsv', sep = '\t', header=TRUE) %>%  
  separate(taxonomy, into = c('Domain', 'Realm', 'Kingdom', 'Phylum', 'Class', 'Order', 'Family'), sep = ';')

contig_abund_pre <- read.table('~/projects/longitudinal_shotgun/data/sqm_tables/19.geNomad_viruses.contigtable', sep = '\t', header = TRUE) 

contig_abund <- contig_abund_pre %>% 
  pivot_longer(values_to = 'TPM', names_to = 'samples', cols = starts_with('TPM')) %>% 
  select(Contig.ID, Tax, TPM, samples) %>% 
  mutate(seq_name = Contig.ID, 
         samples = str_remove_all(samples, 'TPM.')) %>% 
  filter(TPM > 0)


virus <- left_join(viruses, contig_abund, by = 'seq_name') %>% 
  left_join(metadata, by = join_by('samples' == 'Group'))

event_data <- read.table('data/extreme_event_data.csv', sep = ',', header = TRUE)

virus2 <- virus %>% 
  filter(TPM > 0, !is.na(Class), Class != '') %>% 
  group_by(person, day, biota, Class) %>%  
  reframe(sum = sum(TPM)) 

filter(virus2, biota == 'bulk microbiota') %>% 
  ggplot(aes(x = day, y = sum)) +
  geom_rect(data = event_data, aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = Inf, fill = extremevent_type), inherit.aes = FALSE,
            alpha = 0.6) +
  scale_fill_manual(values = c('white','#d94343', '#d98e43', '#f1f011', '#0c9910', '#3472b7', '#7934b7', '#b73485', '#0f5618')) +
  #geom_col() +
  geom_point(aes(color = Class), size =2) +
  geom_line(aes(color = Class), linewidth=1.5) +
  scale_y_log10() +
  facet_wrap(~person, scales = 'free', nrow = 3) +
  labs(x = 'Day', y = '# reads per million [log10]', fill = 'Event')
ggsave('out/geNomad/sumTPM_viruses_bulk.png', dpi = 600)

virus3 <- virus %>% 
  filter(TPM > 0, !is.na(Class), Class != '', Class != 'Caudoviricetes') %>% 
  group_by(person, day, biota, Class) %>%  
  reframe(sum = sum(TPM))

ggplot(virus3, aes(x = day, y = sum)) +
  geom_rect(data = event_data, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = extremevent_type), inherit.aes = FALSE,
            alpha = 0.6) +
  scale_fill_manual(values = c('white','#d94343', '#d98e43', '#f1f011', '#0c9910', '#3472b7', '#7934b7', '#b73485', '#0f5618')) +
  #geom_col() +
  geom_point(aes(color = Class), size =2) +
  geom_line(aes(color = Class), linewidth=1.5) +
  #scale_y_log10() +
  facet_grid(biota~person, scales = 'free') +
  labs(x = 'Day', y = 'Sum TPM', fill = 'Event type')
ggsave('out/geNomad/sumTPM_viruses_allbutCaudoviri.png')


virus4 <- virus %>% 
  filter(TPM > 0, !is.na(Class), Class == 'Caudoviricetes', Family != '') %>% 
  group_by(person, day, biota, Family) %>%  
  reframe(sum = sum(TPM)) %>% 
  filter(biota == 'bulk microbiota')

ggplot(virus4, aes(x = day, y = sum)) +
  geom_rect(data = event_data, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = extremevent_type), inherit.aes = FALSE,
            alpha = 0.6) +
  scale_fill_manual(values = c('white','#d94343', '#d98e43', '#f1f011', '#0c9910', '#3472b7', '#7934b7', '#b73485', '#0f5618')) +
  #geom_col() +
  geom_point(aes(color = Family), size =2) +
  geom_line(aes(color = Family), linewidth=1.5) +
  #scale_y_log10() +
  facet_wrap(~person, scales = 'free', nrow = 3) +
  labs(x = 'Day', y = 'Sum TPM', fill = 'Event type')
ggsave('out/geNomad/sumTPM_viruses_Caudoviricetes.png')


virus_tab <- filter(virus, biota == 'bulk microbiota') %>% 
  select(seq_name, TPM, samples) %>% 
  pivot_wider(names_from = 'seq_name', values_from = 'TPM', values_fill = 0) %>% 
  filter(!is.na(samples)) %>% 
  column_to_rownames('samples')

virus_dist <- vegdist(virus_tab, method = 'bray')
nmds <- metaMDS(virus_dist)                    

nmds_positions <- as.data.frame(scores(nmds, display="sites")) %>%
  rownames_to_column('Group')

dist_meta = nmds_positions %>%
  left_join(metadata, by = 'Group')

dist_meta %>%
  ggplot(aes(x=NMDS1, y=NMDS2, color=person)) +
  geom_point(size = 5) +
  labs(color = 'Individual') +
  theme(legend.position = 'bottom')
ggsave('out/geNomad/nmds.png')

virus_long <- virus_dist %>% 
  as.matrix() %>% 
  as.data.frame() %>% 
  rownames_to_column('Group') %>% 
  pivot_longer(-Group) %>% 
  filter(Group != name) %>% 
  mutate(sample_pairs = paste(Group, name)) %>%
  group_by(sample_pairs) %>%
  reframe(median_value = median(value, na.rm = TRUE)) %>%
  separate(sample_pairs, into = c("Group", "name"), sep = " ") %>%  
  left_join(metadata, by = 'Group') %>% 
  left_join(select(metadata, Group, person, date, biota), by = join_by('name' == 'Group')) %>% 
  mutate(person.x == person.y, 
         same_person = ifelse(person.x == person.y, 'Within individual', 'Between individual'),
         same_biota = ifelse(biota.x == biota.y, 'Same', 'Different'), 
         diff = abs(date.x - date.y)) %>% 
    filter(same_biota == 'Same')

virus_long %>% 
  ggplot(aes(x = same_person, y = median_value, fill = biota.x)) +
  geom_boxplot() +
  scale_fill_manual(values = col) +
  labs(x = '', y = 'Median Bray-Curtis dissimilarity', fill = '') +
  theme(legend.position = 'bottom')
ggsave('out/geNomad/bray_boxplot.png')

virus_long %>% 
  filter(same_person == 'Within individual') %>%  
  ggplot(aes(x = diff, y = median_value, color = biota.x)) +
  geom_point() +
  geom_smooth(se = FALSE) +
  scale_color_manual(values = col) +
  labs(x = 'Days between samples', y = 'Median Bray-Curtis dissimilarity', color = '') +
  theme(legend.position = 'bottom')
ggsave('out/geNomad/bray_time.png')  
