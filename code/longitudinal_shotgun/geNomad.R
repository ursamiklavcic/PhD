# geNomad
library(stringr)
library(ggplot2)
library(readr)
library(dplyr)
library(tibble)
library(lubridate)
library(tidyr)

theme_set(theme_bw())

metadata <- read.table('~/projects/longitudinal_shotgun/data/metadata.csv', header= TRUE, sep = ';') %>%
  mutate(date = dmy(date))

viruses <- read.table('~/projects/longitudinal_shotgun/data/geNomad/01.sqmMicrobiota_virus_summary.tsv', sep = '\t', header=TRUE) %>%  
  separate(taxonomy, into = c('Domain', 'Realm', 'Kingdom', 'Phylum', 'Class', 'Order', 'Family'), sep = ';')

contig_abund_pre <- read.table('~/projects/longitudinal_shotgun/data/sqm_tables/19.geNomad_viruses.contigtable', sep = '\t', header = TRUE) 

contig_abund <- contig_abund_pre %>% 
  pivot_longer(values_to = 'TPM', names_to = 'samples', cols = starts_with('TPM')) %>% 
  select(Contig.ID, Tax, TPM, samples) %>% 
  mutate(seq_name = Contig.ID, 
         samples = str_remove_all(samples, 'TPM.'))


virus <- left_join(viruses, contig_abund, by = 'seq_name') %>% 
  left_join(metadata, by = join_by('samples' == 'Group'))

event_data <- read.table('extreme_event_data.csv', sep = ',', header = TRUE)

virus2 <- virus %>% 
  filter(TPM > 0, !is.na(Class), Class != '') %>% 
  group_by(person, day, biota, Class) %>%  
  reframe(sum = sum(TPM)) 

ggplot(virus2, aes(x = day, y = sum)) +
  geom_rect(data = event_data, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = extremevent_type), inherit.aes = FALSE,
            alpha = 0.6) +
  scale_fill_manual(values = c('white','#d94343', '#d98e43', '#f1f011', '#0c9910', '#3472b7', '#7934b7', '#b73485', '#0f5618')) +
  #geom_col() +
  geom_point(aes(color = Class), size =2) +
  geom_line(aes(color = Class), linewidth=1.5) +
  #scale_y_log10() +
  facet_grid(biota~person, scales = 'free') +
  labs(x = 'Day', y = 'Sum TPM')
ggsave('out/geNomad/sumTPM_viruses.png')

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
  reframe(sum = sum(TPM))

ggplot(virus4, aes(x = day, y = sum)) +
  geom_rect(data = event_data, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = extremevent_type), inherit.aes = FALSE,
            alpha = 0.6) +
  scale_fill_manual(values = c('white','#d94343', '#d98e43', '#f1f011', '#0c9910', '#3472b7', '#7934b7', '#b73485', '#0f5618')) +
  #geom_col() +
  geom_point(aes(color = Family), size =2) +
  geom_line(aes(color = Family), linewidth=1.5) +
  #scale_y_log10() +
  facet_grid(biota~person, scales = 'free') +
  labs(x = 'Day', y = 'Sum TPM', fill = 'Event type')
ggsave('out/geNomad/sumTPM_viruses_Caudoviricetes.png')
