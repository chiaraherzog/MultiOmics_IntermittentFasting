# given dat_wide extracted from multiexperiment, add taxonomic information from count table
# when an ASV is not detected at any timepoint, set value to NA
# author: Charlotte Vavourakis

makeCustomASVtab <- function(dat_wide,sample_type,ASVtab){
  
  require(stringr)
  
  sm_columns = grep(sample_type, names(ASVtab), value = TRUE)
  sm_columns <- sm_columns[!grepl("M0", sm_columns)]
  
  # make a new ASV table from this
  dat_wide = dat_wide %>%
    mutate(primary = paste0(primary,sample_type)) %>%
    column_to_rownames(var="primary") %>%
    t() %>%
    as.data.frame() %>%
    select(sm_columns)
  
  # add taxonomic info
  dat_wide = dat_wide %>%
    filter(rownames(dat_wide) %in% c(ASVtab$OTU))
  taxadd = ASVtab %>% select(OTU, Phylum, Class, Order, Family) %>% dplyr::rename(ASV=OTU)
  sm_columns = grep(sample_type, names(dat_wide), value = TRUE)

  # when an ASV is not detected at any timepoint, set value to NA
  ASV_long = ASVtab %>%
    select(-Kingdom,-Phylum, -Class, -Order, -Family,-Genus,-Species) %>%
    pivot_longer(cols = -c(OTU),
                 names_to = "Sample",
                 values_to = "Count") %>%
    mutate(subjectId = substr(Sample,1,4)) %>%
    mutate(visitId = gsub("SM|FM","",str_sub(Sample,5))) %>%
    group_by(subjectId, OTU) %>%
    mutate(detected = !all(Count == 0)) %>%
    ungroup()
  detected = ASV_long %>%
    select(-subjectId,-visitId,-Count) %>%
    pivot_wider(names_from="Sample", values_from = detected) %>%
    dplyr::rename(ASV=OTU) 
  detected = detected[, c("ASV",sm_columns)] %>%
    column_to_rownames(var="ASV")
  
  if(identical(colnames(detected),colnames(dat_wide))&identical(rownames(detected),rownames(dat_wide))){
    dat_wide[!detected] <- NA
  } else {
    stop("order ASVs in dat_wide, detected don't match. Can't set undetected ASVs to NA")
  }
  
  dat_wide = dat_wide %>%
    rownames_to_column(var="ASV") %>%
    full_join(taxadd)
  dat_wide[sm_columns] = apply(dat_wide[sm_columns], 2, as.numeric)

  return(dat_wide)
  
}