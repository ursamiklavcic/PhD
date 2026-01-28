# EukDetect 

theme_set(theme_bw())

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
