# Determine ethanol resistant species
abund <- read_tsv('~/projects/longitudinal_shotgun/data/metaphlan_abundance_table.txt', comment = '#') %>%
  rename_with(~ str_remove(., '^profiled_'), starts_with('profiled_')) %>%
  #mutate(clade_name = str_remove_all(clade_name, '[a-zA-Z]__')) %>%
  separate(clade_name, into=c('Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species', 'SGB'),
           sep="\\|") %>% 
  mutate(Phylum = ifelse(Phylum == 'p__Firmicutes', 'p__Bacillota', Phylum), 
         Domain = str_remove_all(Domain, 'k__'), 
         Phylum = str_remove_all(Phylum, 'p__'), 
         Class = str_remove_all(Class, 'c__'), 
         Order = str_remove_all(Order, 'o__'), 
         Family = str_remove_all(Family, 'f__'), 
         Genus = str_remove_all(Genus, 'g__'), 
         Species = str_remove_all(Species, 's__'), 
         SGB = str_remove_all(SGB, 't__')) %>% 
  select(-MC013)

bacteria <- filter(abund, Domain == 'Bacteria', !is.na(Phylum), !is.na(Class), 
                   !is.na(Order), !is.na(Family), !is.na(Genus), !is.na(Species), !is.na(SGB)) %>% 
  pivot_longer(-c(Domain, Phylum, Class, Order, Family, Genus, Species, SGB)) %>% 
  left_join(metadata, by = join_by('name' == 'Group')) %>% 
  mutate(PA = ifelse(value > 0, 1, 0)) %>% 
  select(Domain, Phylum, Class, Order, Family, Genus, Species, SGB, original_sample, value, name)

etoh_species <- full_join(bacteria %>% filter(substr(name, 1, 1) == 'M'), 
                          bacteria %>% filter(substr(name, 1, 1) == 'S'), 
                          by = join_by('Domain', 'Phylum', 'Class', 'Order', 'Family', 
                                       'Genus', 'Species', 'SGB', 'original_sample')) %>%
  # Define if species in a sample of stool is ethanol resistant 
  # Contition 1: present in both bulk microbiota sample and ethanol resistant fraction
  # Condition 2: higher relative abudnance in EtOH sample than microbiota
  mutate(is_etoh_resistant = ifelse(value.x > 0 & value.y > 0 & value.y > value.x, 'Yes', 'No')) %>%
  group_by(SGB) %>%
  # Calculate the number of times a species was present in samples
  reframe(no_present = n_distinct(name.y, na.rm = TRUE), 
          # Caluclate how many times OTU was defined as part of EtOH fraction based on Conditions 1 & 2
          no_Yes = ceiling(sum(is_etoh_resistant == 'Yes', na.rm = TRUE))) %>%
  # Species that have been defined as part of the ethanol resistant fraction in at least 5% of samples where they were found! 
  # (to avoid mistakes of protocol and exclude highly abundant species that maybe were seen as ethanol resistant but just didn't get destoryed!)
  filter(no_Yes > (no_present * 0.05)) %>%
  pull(unique(SGB))

species_filt <- filter(bacteria, SGB %in% etoh_species) %>% 
  select(Species, SGB) %>% 
  distinct() %>% 
  mutate(is_etoh_resistant = TRUE)
write_tsv(species_filt, 'data/longitudinal_shotgun/ethanol_resistant_SGB.tsv')
