# EukDetect 

theme_set(theme_bw(base_size=14))

hits <- read.table('../longitudinal_shotgun/data/EukDetect/all_hits.txt', sep = '\t', header = T, fill=TRUE, quote="")

read_counts <- read.table('../longitudinal_shotgun/data/EukDetect/read_counts.txt', sep = '\t', header = F) %>% 
  rename('sample' = 'V1', 'reads' = 'V2') %>% 
  mutate(sample = substr(sample, 1, 5)) %>%  
  group_by(sample) %>% 
  reframe(reads = mean(reads))

hits %>% 
  group_by(present) %>% 
  reframe(n = n_distinct(sample))

# # A tibble: 3 × 2
# present     n
# <chr>   <int>
# 1 YES        68
# 2 pass       13
# 3 NA         35

hits_yes <- hits %>% filter(present == 'YES')

hits_yes %>%  
  left_join(read_counts, 'sample') %>% 
  mutate(percent_reads = Read_counts/reads*100, 
         person = substr(sample, 2, 2)) %>% 
  ggplot(aes(x = percent_reads, y = sample, fill = Name)) +
  geom_col() +
  labs(x = 'Percentage of all reads', y = '', fill= 'Organism') +
  facet_wrap(~person, scales = 'free')

ggsave('out/EukDetect/species_percent_reads.svg')

# A bit different plot 
metadata <- read.table('data/metadata.csv', sep = ';', header = T)

event_data <- read.table('data/extreme_event_data.csv', sep = ',', header = TRUE) %>%  
  filter(person %in% c('B', 'C', 'D', 'E', 'F', 'G', 'H', 'I'))


hits_yes %>%  
  left_join(read_counts, 'sample') %>% 
  left_join(metadata, by = join_by('sample' == 'Group')) %>% 
  mutate(percent_reads = Read_counts/reads*100, 
         biota = ifelse(biota == 'bulk microbiota', 'untreated sample', biota), 
         Name = ifelse(grepl('Blastocystis', Name), 'Blastocystis sp.', Name)) %>% 
  ggplot(aes(x = day, y = percent_reads)) +
  geom_rect(data = event_data, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = extremevent_type), inherit.aes = FALSE,
            alpha = 0.6) +
  scale_fill_manual(values = c('white','#d94343', '#d98e43', '#f1f011', '#0c9910', '#3472b7', '#7934b7', '#b73485', '#0f5618')) +
  geom_point(size = 5, aes(color = Name)) +
  geom_line(linewidth = 2, aes(color = Name)) +
  #scale_y_continuous(limits = c(0, NA), expand = c(0, 0)) +
  labs(y = 'Percentage of all reads [%]', x = 'Day', fill= 'Event', color = 'Organism') +
  facet_wrap(~person, scales = 'free')

ggsave('out/EukDetect/plot.svg')

# Abundant species 
hits_yes %>%  
  left_join(read_counts, 'sample') %>% 
  mutate(percent_reads = Read_counts/reads*100, 
         person = substr(sample, 2, 2)) %>% 
  group_by(Name) %>% 
  reframe(sum = sum(percent_reads))

hits_yes %>% left_join(read_counts, 'sample') %>% 
  group_by(present) %>% 
  reframe(all = sum(Read_counts)/sum(reads)*100)
