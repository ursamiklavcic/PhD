# Properties of OTUs found in only ethanol treated samples 

library(readr)
library(tidyr)
library(dplyr)
library(tibble)
library(vegan)
library(lubridate)
library(ggpubr)
library(scales)

set.seed(96)
theme_set(theme_bw(base_size=14))

otutab <- readRDS('data/longitudinal_amplicons/otutab_ethanol_bulk.RDS')
taxtab <- readRDS('data/longitudinal_amplicons/taxtab.RDS')
metadata <- readRDS('data/longitudinal_amplicons/metadata.RDS')
ddPCR <- readRDS('data/longitudinal_amplicons/ddPCR.RDS')


otu_long <- rownames_to_column(as.data.frame(otutab), 'Group') %>% 
  pivot_longer(cols = starts_with('Otu')) %>%
  left_join(metadata %>% select(original_sample, Group, person, date), by = 'Group') %>%
  group_by(Group) %>%
  mutate(rel_abund = value / sum(value), 
         PA = ifelse(value > 0, 1, 0)) %>%
  ungroup() %>%
  left_join(ddPCR, by = join_by('Group' == 'Sample')) %>%
  mutate(norm_abund = rel_abund * copies) %>%
  select(Group, name, value, original_sample, person, norm_abund, rel_abund, PA, date) %>%
  left_join(taxtab, by = 'name')

# Filter OTUs only in EtOH samples 
# OTUs found only in ethanol treated samples 

otus_in_both <- otu_long %>%
  filter(PA == 1) %>%  
  pull(unique(name))

otus_only_untreated <- otu_long %>%
  filter(substr(Group, 1, 1) == 'M' & PA == 1) %>% 
  pull(unique(name))

otus_etoh_only <- otu_long %>% filter(name %in% otus_in_both) %>%  
  filter(!name %in% otus_only_untreated) %>% 
  filter(substr(Group, 1, 1) == 'S' & value > 0) 

otus_etoh_only %>% reframe(n = n_distinct(name))

otus_etoh_only %>% 
  ggplot(aes(x = rel_abund, y = Phylum, fill = Phylum)) +
  geom_boxplot() +
  scale_x_log10() +
  labs(x = 'Relative abundance [log10]', y = '')
ggsave('out/only_etoh_OTUs/relative_abundnace_phylum.png', dpi=600)

# Number per phylum 
otus_etoh_only %>% 
   group_by(Phylum) %>%  
   reframe(n = n_distinct(name)) %>% 
   ggplot(aes(x = n, y = Phylum, fill = Phylum)) +
   geom_col() +
   geom_label(aes(label = n)) +
  labs(x = '# OTUs', y = '')
