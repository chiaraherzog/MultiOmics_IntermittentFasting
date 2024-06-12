#' @name bubble_summary
#' @description
#' Creates df for bubble plot
#' @param lmmdat LMM data frame (minimal model)
#' @param t_test T test (overall)
#' @param vars list of variables to be included
#' @return Returns a data frame for the bubble diagram

bubble_summary <- function(lmmdat, t_test, vars){
  
  # Linear mixed models: main data, excluding metabolome
  dat <- lmmdat |> 
    tidyr::separate(x, into = c("assay", "Feature"), sep = "_", extra = "merge") |> 
    dplyr::inner_join(vars, by = c('assay' = 'assay',
                                   'Feature' = 'x'))|> 
    dplyr::mutate(assay = case_when(((assay == "Blood haemogram") | assay2 == 'Blood test') ~ 'Routine bloods',
                                    (grepl('methylation', assay) & grepl('cervical', assay)) ~ 'Cervical methylation',
                                    (grepl('methylation', assay) & grepl('buccal', assay)) ~ 'Buccal methylation',
                                    (grepl('methylation', assay) & grepl('blood', assay)) ~ 'Blood methylation',
                                    grepl("Skin", assay) ~ "Skin histology and TEWL",
                                    grepl("vifat|scfat", Feature) ~ "Subcutaneous and visceral fat",
                                    !is.na(assay2) ~ assay2,
                                    TRUE ~ assay)) |> 
    dplyr::select(assay, Feature, contains("p.value_visitId")) |> 
    dplyr::mutate(assay = gsub(": families", "", assay)) |> 
    dplyr::filter(!grepl("ASVs|metabolome", assay))
  
  dat_m <- t_test |> 
    tidyr::separate(rowname, into = c("assay", "Feature"), sep = "_", extra = "merge") |> 
    dplyr::inner_join(vars, by = c('assay' = 'assay',
                                   'Feature' = 'x'))|> 
    dplyr::filter(grepl("nuclear", assay)) |> 
    dplyr::mutate(assay = assay2) |> 
    dplyr::select(assay, Feature, starts_with("p_")) |> 
    dplyr::rename_with(~ gsub("p_", "p.value_visitId", .x), starts_with("p_"))
  
  dat <- dat |> 
    dplyr::bind_rows(dat_m) |> 
    pivot_longer(starts_with("p.value_visitIdM"), names_to = "visitId", values_to = "p.value") |> 
    mutate(significant = ifelse(p.value < 0.05, 1, 0)) |> 
    group_by(assay, visitId) |> 
    summarise(
      num_significant = sum(significant),
      proportion_significant = sum(significant) / n()) |> 
    ungroup()
  dat$visitId <- gsub("p.value_visitId","",dat$visitId)
  
  return(dat)
}