# given an ASV count table and a long format extracted from multi-assay experiment, relabel taxon level of interest
# in for example /data/data_baseline_change.Rdata there are some dots replacing spaces and -, so this script won't work here!
# e.g. Family; ASV or Phylum; Family
# Author: Charlotte Vavourakis

relabel_microbiomefeat <- function(dat_long,datASV,level,leveladd){
  
  relab <- datASV %>%
    dplyr::rename(ASV = OTU) %>%
    mutate(label = paste(get(leveladd),get(level),sep = ";")) %>%
    select(level,label) %>%
    distinct()
  
  relab <- dat_long %>%
    full_join(relab) %>%
    select(-level) %>%
    select(label, everything())
  
  return(relab)
  
}