#' @name loadRMcorr
#' @description
#' Loads and prepares rmcorr data
#' @param filter_ASV Should ASVs be removed? False by default.
#' @return returns correlation matrix.

loadRMcorr <- function(filter_ASV = F){
  
  suppressPackageStartupMessages(library(dplyr))
  # Loading repeated measures correlation - raw and change values
  load(here('out/rmcorr_change.Rdata'))
  rmcorr_change <- rmcorr
  load(here('out/rmcorr.Rdata'))
  
  # First: filter those that have consistent sign for change
  rm <- rmcorr |> 
    dplyr::select(measure1, measure2, assay1, assay2, rmcorr.r, p.vals) |> 
    dplyr::left_join(dplyr::select(rmcorr_change, measure1, measure2, rmcorr.r, p.vals, assay1, assay2),
                     by = c('measure1', 'measure2', 'assay1', 'assay2'), suffix = c('', '.change')) |> 
    dplyr::filter(sign(rmcorr.r) == sign(rmcorr.r.change) | grepl("-log", assay1) | grepl("-log", assay2)) |> 
    
    # Additional filtering: Exclude visitId
    dplyr::filter(!assay1 %in% c('visitId') & !assay2 %in% c('visitId')) |>
    dplyr::filter(!grepl("Cotinine", assay1) & !grepl("Cotinine", assay2))
  
  # Relabel assays
  assaylabels <- c("Flow cytometry: white blood cell staining" = "Flow cytometry:\nimmune cells",
                   "Flow cytometry: T cell staining" = "Flow cytometry:\nimmune cells",
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
    dplyr::mutate(assay1 = recode(assay1, !!! assaylabels),
                  assay2 = recode(assay2, !!! assaylabels))
  corr <- corr |> 
    dplyr::mutate(assay1 = case_when((grepl("exam", assay1) & !grepl("_fe|_fv|_sysbp|_diabp|_vo2max|_rel|_abs|_max", measure1)) ~ "Routine\nbloods",
                                     ((grepl("exam", assay1) & grepl("_fe|_fv|_sysbp|_diabp|_vo2max|_rel|_abs|_max", measure1)) | grepl("pwv|imt|plaque", measure1)) ~ "Functional\nclinical measures",
                                     grepl("bmi|weight|scfat|vifat|bcm|ecw|fm", measure1) ~ "Body\ncomposition",
                                     TRUE ~ assay1),
                  assay2 = case_when((grepl("exam", assay2) & !grepl("_fe|_fv|_sysbp|_diabp|_vo2max|_rel|_abs|_max", measure2)) ~ "Routine\nbloods",
                                     ((grepl("exam", assay2) & grepl("_fe|_fv|_sysbp|_diabp|_vo2max|_rel|_abs|_max", measure2)) | grepl("pwv|imt|plaque", measure2)) ~ "Functional\nclinical measures",
                                     grepl("bmi|weight|scfat|vifat|bcm|ecw|fm", measure2) ~ "Body\ncomposition",
                                     TRUE ~ assay2))
  
  # Filter obvious 'self correlations'
  corr <- corr |>
    dplyr::filter(assay1 != assay2 & !(assay1 == 'Immune\nage' & grepl("cytometry", assay2)) & !(assay2 == 'Immune\nage' & grepl("cytometry", assay1))) |>
    dplyr::filter(assay1 != assay2 & !(assay1 == 'ImmuneSMK' & grepl("cytometry", assay2)) & !(assay2 == 'ImmuneSMK' & grepl("cytometry", assay1))) |>
    dplyr::filter(!(grepl("Saliva\nmicrobiome", assay1) & grepl("Saliva\nmicrobiome", assay2))) |> 
    dplyr::filter(!(grepl("Stool\nmicrobiome", assay1) & grepl("Stool\nmicrobiome", assay2))) |> 
    dplyr::filter(!(grepl("Saliva\nmicrobiome", assay1) & grepl("Saliva microbiome: ASVs", assay2))) |> 
    dplyr::filter(!(grepl("Stool\nmicrobiome", assay1) & grepl("Stool microbiome: ASVs", assay2))) |> 
    dplyr::filter(!(grepl("Saliva\nmetabolome", assay1) & grepl("Saliva\nmetabolome", assay2))) |> 
    dplyr::filter(!(grepl("Urine\nmetabolome", assay1) & grepl("Urine\nmetabolome", assay2)))
  
  # Filter p value
  corr <- corr |> 
    # dplyr::filter(p.vals < 0.01) |> 
    dplyr::filter(!(grepl("metabolome", assay1) & !grepl("log", assay1)) & !(grepl("metabolome", assay2) & !grepl("log", assay2)))
  
  
  if(filter_ASV==T){
    corr <- corr |> 
      dplyr::filter(!grepl("ASV", assay1) & !grepl("ASV", assay2))
  }
  
  corr <- corr |> 
    dplyr::filter(!is.na(p.vals)) |> 
    dplyr::arrange(p.vals) |> 
    dplyr::mutate(padj = p.adjust(p.vals, method = 'fdr'))
  
  return(corr)
}
