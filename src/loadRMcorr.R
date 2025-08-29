#' @name loadRMcorr
#' @description
#' Loads and prepares rmcorr data
#' @param filter_ASV Should ASVs be removed? False by default.
#' @return returns correlation matrix.

loadRMcorr <- function(filter_ASV = F){
  
  suppressPackageStartupMessages(library(dplyr))
  # Loading repeated measures correlation - raw and change values
  load(here('out/rmcorr_change.Rdata'))
  rmcorr_change <- res_change
  load(here('out/rmcorr.Rdata'))
  rmcorr <- res
  
  # First: filter those that have consistent sign for change
  rm1 <- rmcorr |> 
    dplyr::filter(measure1 != measure2) |> 
    tidyr::separate(measure1, sep = "_", into = c("assay1", "measure1"), remove = T,
                    extra = 'merge') |> 
    tidyr::separate(measure2, sep = "_", into = c("assay2", "measure2"), remove = T,
                    extra = 'merge')
  
  rm2 <- rmcorr_change |> 
    dplyr::filter(measure1 != measure2) |> 
    tidyr::separate(measure1, sep = "_", into = c("assay1", "measure1"), remove = T,
                    extra = 'merge') |> 
    tidyr::separate(measure2, sep = "_", into = c("assay2", "measure2"), remove = T,
                    extra = 'merge')
  
  rm <- rm1 |> 
    dplyr::select(measure1, measure2, assay1, assay2, r, p) |> 
    dplyr::left_join(dplyr::select(rm2, measure1, measure2, r, p, assay1, assay2),
                     by = c('measure1', 'measure2', 'assay1', 'assay2'), suffix = c('', '.change')) |> 
    dplyr::filter(sign(r) == sign(r.change)) |> 
    
    # Additional filtering: Exclude visitId
    dplyr::filter(!assay1 %in% c('visitId') & !assay2 %in% c('visitId'))
  
  # Relabel assays
  assaylabels <- c("Flow cytometry: white blood cell staining" = "Flow cytometry:\nimmune cells",
                   "Flow cytometry: T cell staining" = "Flow cytometry:\nimmune cells",
                   "Flow cytometry: stimulated T cells" = "Flow cytometry:\nimmune cell stimulation",
                   "Flow cytometry: unstimulated T cells" = "Flow cytometry:\nimmune cell stimulation",
                   "Blood haemogram" = "Routine\nbloods",
                   "Skin histology and transepidermal water loss assay" = "Functional\nclinical measures",
                   "Body composition" = "Body\ncomposition",
                   "Vascular and body sonography" = "Functional\nclinical measures",
                   "Composite methylation scores: buccal" = "Buccal\nmethylation",
                   "Composite methylation scores: cervical" = "Cervical\nmethylation",
                   "Composite methylation scores: blood" = "Blood\nmethylation",
                   "Urine nuclear magnetic resonance: normalized" = "Urine\nmetabolome",
                   "Saliva nuclear magnetic resonance: normalized" = "Saliva\nmetabolome",
                   "Saliva nuclear magnetic resonance: normalized-log" = "Saliva\nmetabolome-log",
                   "Urine nuclear magnetic resonance: normalized-log" = "Urine\nmetabolome-log",
                   "Saliva microbiome: families" = "Saliva\nmicrobiome",
                   "Stool microbiome: families" = "Stool\nmicrobiome",
                   "Immune age: general" = "Immune\nage",
                   "ImmuneSMK" = 'ImmuneSMK')
  
  
  corr <- rm |> 
    dplyr::mutate(assay1_relab = recode(assay1, !!! assaylabels),
                  assay2_relab = recode(assay2, !!! assaylabels))
  
  corr <- corr |> 
    dplyr::mutate(assay1_relab = case_when((grepl("exam", assay1_relab) & !grepl("_fe|_fv|_sysbp|_diabp|_vo2max|_rel|_abs|_max", measure1)) ~ "Routine\nbloods",
                                     ((grepl("exam", assay1_relab) & grepl("_fe|_fv|_sysbp|_diabp|_vo2max|_rel|_abs|_max", measure1)) | grepl("pwv|imt|plaque", measure1)) ~ "Functional\nclinical measures",
                                     grepl("bmi|weight|scfat|vifat|bcm|ecw|fm", measure1) ~ "Body\ncomposition",
                                     TRUE ~ assay1_relab),
                  assay2_relab = case_when((grepl("exam", assay2_relab) & !grepl("_fe|_fv|_sysbp|_diabp|_vo2max|_rel|_abs|_max", measure2)) ~ "Routine\nbloods",
                                     ((grepl("exam", assay2_relab) & grepl("_fe|_fv|_sysbp|_diabp|_vo2max|_rel|_abs|_max", measure2)) | grepl("pwv|imt|plaque", measure2)) ~ "Functional\nclinical measures",
                                     grepl("bmi|weight|scfat|vifat|bcm|ecw|fm", measure2) ~ "Body\ncomposition",
                                     TRUE ~ assay2_relab))
  
  # Filter obvious 'self correlations'
  corr <- corr |>
    dplyr::filter(assay1_relab != assay2_relab & !(assay1_relab == 'Immune\nage' & grepl("cytometry", assay2_relab)) & !(assay2_relab == 'Immune\nage' & grepl("cytometry", assay1_relab))) |>
    dplyr::filter(assay1_relab != assay2_relab & !(assay1_relab == 'ImmuneSMK' & grepl("cytometry", assay2_relab)) & !(assay2_relab == 'ImmuneSMK' & grepl("cytometry", assay1_relab))) |>
    dplyr::filter(!(grepl("Saliva\nmicrobiome", assay1_relab) & grepl("Saliva\nmicrobiome", assay2_relab))) |> 
    dplyr::filter(!(grepl("Stool\nmicrobiome", assay1_relab) & grepl("Stool\nmicrobiome", assay2_relab))) |> 
    dplyr::filter(!(grepl("Saliva\nmicrobiome", assay1_relab) & grepl("Saliva microbiome: ASVs", assay2_relab))) |> 
    dplyr::filter(!(grepl("Stool\nmicrobiome", assay1_relab) & grepl("Stool microbiome: ASVs", assay2_relab))) |> 
    dplyr::filter(!(grepl("Saliva\nmetabolome", assay1_relab) & grepl("Saliva\nmetabolome", assay2_relab))) |> 
    dplyr::filter(!(grepl("Urine\nmetabolome", assay1_relab) & grepl("Urine\nmetabolome", assay2_relab)))
  
  # Filter p value
  corr <- corr |> 
    # dplyr::filter(p.vals < 0.01) |> 
    dplyr::filter(!(grepl("metabolome", assay1_relab) & !grepl("log", assay1_relab)) & !(grepl("metabolome", assay2_relab) & !grepl("log", assay2_relab)))
  
  
  if(filter_ASV==T){
    corr <- corr |> 
      dplyr::filter(!grepl("ASV", assay1_relab) & !grepl("ASV", assay2_relab))
  }
  
  # Remove duplicates
  corr <- corr |> 
    dplyr::rowwise() |> 
    mutate(
      # Create canonical ordering: paste the "smaller" first
      pair1 = paste(assay1, measure1, sep = "|"),
      pair2 = paste(assay2, measure2, sep = "|"),
      key   = paste(sort(c(pair1, pair2)), collapse = "___")  # order-independent ID
    ) |> 
    ungroup() |> 
    distinct(key, .keep_all = TRUE) |> 
    dplyr::select(-pair1, -pair2, -key) 
  
  corr <- corr |> 
    dplyr::mutate(new = paste0(assay1_relab, measure1)) |> 
    dplyr::filter(!is.na(p)) |> 
    dplyr::arrange(p) |> 
    dplyr::mutate(padj = p.adjust(p, method = 'fdr')) |> 
    dplyr::group_by(new) |> 
    dplyr::mutate(padj_v2 = p.adjust(p, method = 'fdr')) |>
    dplyr::ungroup() |> 
    dplyr::select(-new)
  
  # Remove
  
  return(corr)
}
