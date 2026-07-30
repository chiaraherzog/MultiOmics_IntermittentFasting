library(ampvis2)
library(tidyverse)

mean_rel_abund_visit <- function(ASVtable, pheno, agg = "Family") {
  
  # load and normalise to relative abundance
  
  d <- amp_load(otutable = ASVtable, metadata = pheno)
  d <- normaliseTo100(d)
  
  # aggregate to family level
  
  d <- aggregate_abund(
    d$abund,
    d$tax,
    tax_aggregate = agg,
    format = "abund",
    calcSums = FALSE
  )
  
  # convert to long format and join metadata
  
  d_long <- d %>%
    t() %>%
    as.data.frame() %>%
    tibble::rownames_to_column("sampleId") %>%
    pivot_longer(
      cols = -sampleId,
      names_to = "tax",
      values_to = "rel_abundance"
    ) %>%
    left_join(pheno %>% select(sampleId, visitId), by = "sampleId")
  
  # compute mean relative abundance per family per visit
  
  summary <- d_long %>%
    group_by(tax, visitId) %>%
    summarise(
      mean_rel_abundance = mean(rel_abundance, na.rm = TRUE),
      sd_rel_abundance = sd(rel_abundance, na.rm = TRUE),
      .groups = "drop"
    )
  
  return(summary)
}
