# File of bacterial species for sporulation determination
spore_determine <-  read_tsv('~/projects/longitudinal_shotgun/data/metaphlan_abundance_table.txt', comment = '#') %>%
  rename_with(~ str_remove(., '^profiled_'), starts_with('profiled_')) %>% 
  filter(grepl('s__', clade_name), grepl('Bacteria', clade_name), !grepl('t__', clade_name)) %>% 
  select(clade_name) %>%  
  mutate(Species = sapply(strsplit(clade_name, "\\|"), function(x) {
    sp <- x[grep("^s__", x)]
    if (length(sp) == 0) return(NA_character_)
    sub("^s__", "", sp)}), 
    Species = str_replace_all(Species, '_', ' '))

write_tsv(spore_determine, 'data/sporulation_ability/spore_determine.tsv')

# Find reference genomes for all species! - script find_genomes_from_name.sh with datasets conda env

