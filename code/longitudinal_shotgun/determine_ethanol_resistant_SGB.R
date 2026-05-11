# Determine ethanol resistant species

library(microbiomics)

abund <- read_metaphlan_table('~/projects/longitudinal_shotgun/data/metaphlan_abundance_table.txt', kingdom = "k__Bacteria", 
                              lvl = 7, normalize = TRUE) %>% 
  rownames_to_column('name') %>% 
  pivot_longer(names_to = 'tax', values_to = 'value', cols = starts_with('k__')) %>% 
  mutate(name = str_remove_all(name, 'profiled_')) %>% 
  tidyr::separate_wider_delim(tax, delim = ".",
                              names = c('Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species')) %>% 
  mutate(Phylum = ifelse(Phylum == 'p__Firmicutes', 'p__Bacillota', Phylum), 
         Domain = str_remove_all(Domain, 'k__'), 
         Phylum = str_remove_all(Phylum, 'p__'), 
         Class = str_remove_all(Class, 'c__'), 
         Order = str_remove_all(Order, 'o__'), 
         Family = str_remove_all(Family, 'f__'), 
         Genus = str_remove_all(Genus, 'g__'), 
         Species = str_remove_all(Species, 's__')) %>% 
  filter(name != 'MC013') %>% 
  left_join(metadata, by = join_by('name' == 'Group')) %>% 
  mutate(PA = ifelse(value > 0, 1, 0)) %>% 
  select(Domain, Phylum, Class, Order, Family, Genus, Species, original_sample, value, name)

etoh_species <- full_join(abund %>% filter(substr(name, 1, 1) == 'M'), 
                          abund %>% filter(substr(name, 1, 1) == 'S'), 
                          by = join_by('Domain', 'Phylum', 'Class', 'Order', 'Family', 
                                       'Genus', 'Species', 'original_sample')) %>%
  # Define if species in a sample of stool is ethanol resistant 
  # Contition 1: present in both bulk microbiota sample and ethanol resistant fraction
  # Condition 2: higher relative abudnance in EtOH sample than microbiota
  mutate(is_etoh_resistant = ifelse(value.x > 0 & value.y > 0 & value.y > value.x, 'Yes', 'No')) %>%
  group_by(Species) %>%
  # Calculate the number of times a species was present in samples
  reframe(no_present = n_distinct(name.y, na.rm = TRUE), 
          # Caluclate how many times OTU was defined as part of EtOH fraction based on Conditions 1 & 2
          no_Yes = ceiling(sum(is_etoh_resistant == 'Yes', na.rm = TRUE))) %>%
  # Species that have been defined as part of the ethanol resistant fraction in at least 5% of samples where they were found! 
  # (to avoid mistakes of protocol and exclude highly abundant species that maybe were seen as ethanol resistant but just didn't get destoryed!)
  filter(no_Yes > (no_present * 0.05)) %>%
  pull(unique(Species))
length(unique(etoh_species))

only_etoh_species <- full_join(abund %>% filter(substr(name, 1, 1) == 'M'),
                               abund %>% filter(substr(name, 1, 1) == 'S'), 
                               by = join_by('Domain', 'Phylum', 'Class', 'Order', 'Family', 
                                            'Genus', 'Species','original_sample')) %>%
  filter(value.x == 0, !is.na(value.y), value.y > 0) %>% 
  pull(unique(Species))
length(unique(only_etoh_species))


uncertain_species <- full_join(abund %>% filter(substr(name, 1, 1) == 'M'), 
                               abund %>% filter(substr(name, 1, 1) == 'S'), 
                               by = join_by('Domain', 'Phylum', 'Class', 'Order', 'Family', 
                                            'Genus', 'Species', 'original_sample')) %>%
  mutate(is_etoh_resistant = ifelse(value.x > 0 & value.y > 0 & value.y > value.x, 'Yes', 'No')) %>%
  group_by(Species) %>%
  reframe(no_present = n_distinct(name.y, na.rm = TRUE), 
          no_Yes = ceiling(sum(is_etoh_resistant == 'Yes', na.rm = TRUE))) %>%
  # Filter OTUs that were detected as EtOH resistant at least once, but were detected as such in less than 5% of samples, to exclude them from the analysis 
  filter(no_Yes > 1) %>%
  filter(no_Yes < (no_present * 0.05)) %>%
  filter(!Species %in% only_etoh_species) %>% 
  pull(unique(Species))
length(unique(uncertain_species))


etoh_species_table <- abund %>% 
  mutate(is_ethanol_resistant = case_when(
    Species %in% etoh_species ~ "Ethanol-resistant",
    Species %in% only_etoh_species ~ "Only ethanol treated samples",
    Species %in% uncertain_species ~ "Uncertain",
    TRUE ~ "Non ethanol-resistant" )) %>% 
  select(Domain, Phylum, Class, Order, Family, Genus, Species, is_ethanol_resistant) %>% 
  distinct()

write_tsv(etoh_species_table, '~/projects/thesis/data/longitudinal_shotgun/ethanol_resistant_species.tsv')
  
length(unique(filter(etoh_species_table, is_ethanol_resistant == "Non ethanol-resistant")$Species)) 

  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
