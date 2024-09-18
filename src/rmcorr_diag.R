#' @name multi_network_diagram
#' @param features Features that should be included in the diagram
#' @param seed set for the layout
#' @param layout layout option
#' @param legend should a lgend be printed?

rmcorr_diag <- function(feature){
  
  # extrace corrs with that feature
  tmp <- corr |> 
    dplyr::filter(grepl(feature, measure1) | grepl(feature, measure2)) |>
    dplyr::mutate(measure = ifelse(grepl(feature, measure1), measure2, measure1),
                  assay = stringr::str_split_i(measure, "_", 1)) |> 
    dplyr::rowwise() |> 
    dplyr::mutate(feature = gsub(paste0(assay, "_"), "", measure)) |> 
    dplyr::ungroup()
  
  tmp$label <-  vars[match(tmp$feature, vars$x),]$label
  tmp <- tmp |>
    dplyr::mutate(label = case_when(is.na(label) & grepl("microbiome", assay) ~ feature,
                                    is.na(label) & grepl("nuclear", assay) ~ gsub("[.]", " ", gsub("^X", "", feature)),
                                     TRUE ~ label))
  
  # Relabel assays Assay labels for the diagram
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
                   "Saliva nuclear magnetic resonance: normalized-log" = "Saliva\nmetabolome",
                   "Urine nuclear magnetic resonance: normalized-log" = "Urine\nmetabolome",
                   "Saliva microbiome: families" = "Saliva\nmicrobiome",
                   "Stool microbiome: families" = "Stool\nmicrobiome",
                   "Immune age: general" = "Immune\nage")
  
  tmp <- tmp |> 
    dplyr::mutate(assay = recode(assay, !!! assaylabels)) |> 
    dplyr::mutate(assay = case_when(grepl("cholesterol|bilirubin|transferrin|urea|ldh|ferritin|uricacid|vitb9|creatinine|hdl", feature) ~ "Routine\nbloods",
                                    grepl("systolic|VO2|FEV|FVC|power", label) ~ "Functional\nclinical measures",
                                    TRUE ~ assay))
  
  featurelabel = vars[grepl(feature, vars$x),]$label
  
  p <- tmp |> 
    ggplot(aes(x = rmcorr.r,
               y = forcats::fct_reorder(label, rmcorr.r),
               size = -log10(p.vals),
               colour = assay)) +
    geom_point(alpha = 0.9) +
    geom_vline(xintercept = 0, linetype = 'dotted') +
    coord_cartesian(xlim = c(-1, 1)) +
    theme_bw() +
    theme(axis.title.x = element_markdown(),
          axis.ticks.y = element_blank()) +
    scale_colour_manual(values = grid.col,
                        name = '') +
    scale_size_continuous(range = c(2, 6),
                          limits = c(2, 25),
                          name = '-log10(p value)') +
    labs(x = 'r<sub>rm</sub>',
         y = '',
         subtitle = featurelabel)
  
  return(p)
  
}
