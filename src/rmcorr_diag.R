#' @name multi_network_diagram
#' @param features Features that should be included in the diagram
#' @param seed set for the layout
#' @param layout layout option
#' @param legend should a lgend be printed?

rmcorr_diag <- function(feature, corrObj){
  
  # extrace corrs with that feature
  tmp <- corrObj |> 
    dplyr::filter(grepl(feature, measure1) | grepl(feature, measure2)) |>
    dplyr::mutate(measure_cor = ifelse(grepl(feature, measure1), label2, label1),
                  assay_cor = ifelse(grepl(feature, measure1), assay2, assay1)) 
  
  # Relabel assays for this specific plot
  assaylabels <- c("Flow cytometry: white blood cell staining" = "Flow cytometry:\nimmune cells",
                   "Flow cytometry: T cell staining" = "Flow cytometry:\nimmune cells",
                   "Flow cytometry: stimulated T cells" = "Flow cytometry:\nimmune cell stimulation",
                   "Flow cytometry: unstimulated T cells" = "Flow cytometry:\nimmune cell stimulation",
                   "Blood haemogram" = "Routine\nbloods",
                   "Skin histology and transepidermal water loss assay" = "Clinical\nmeasures",
                   "Body composition" = "Body\ncomposition",
                   "Vascular and body sonography" = "Clinical\nmeasures",
                   'Functional sports exam' = "Clinical\nmeasures",
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
  
  col12 <- grDevices::colorRampPalette(cols[c(8, 1,2,3,5,4,6,7)])(13)
  grid.col = c("Flow cytometry:\nimmune cells" = col12[11],
               "Immune\nage" = col12[10],
               "Routine\nbloods" = col12[9],
               "Body\ncomposition" = col12[8],
               "Clinical\nmeasures" = col12[7],
               "Cervical\nmethylation" = col12[6],
               "Buccal\nmethylation" = col12[5],
               "Blood\nmethylation" = col12[4],
               "Saliva\nmicrobiome" = col12[3],
               "Stool\nmicrobiome" = col12[2],
               "Urine\nmetabolome" = col12[1],
               "Saliva\nmetabolome" = col12[12],
               "Flow cytometry:\nimmune cell stimulation" = col12[13])
  
  tmp <- tmp |> 
    dplyr::mutate(assay_cor = recode(assay_cor, !!! assaylabels))
  
  # fix immune labels 
  tmp <- tmp |> 
    dplyr::mutate(measure_cor = gsub("[-]", "<sup>-</sup>", measure_cor),
                  measure_cor = gsub("[+]", "<sup>+</sup>", measure_cor),
                  measure_cor = gsub("terminally differentiated", "term. diff.", measure_cor),
                  measure_cor = gsub("VO2peak", "VO<sub>2</sub>peak", measure_cor)
    )
  
  p <- tmp |> 
    ggplot(aes(x = r,
               y = forcats::fct_reorder(measure_cor, r),
               size = -log10(p),
               colour = assay_cor)) +
    geom_point(alpha = 0.9) +
    geom_vline(xintercept = 0, linetype = 'dotted') +
    coord_cartesian(xlim = c(-1, 1)) +
    theme_bw() +
    theme(axis.title.x = element_markdown(),
          axis.ticks.y = element_blank(),
          axis.text.y = element_markdown()) +
    scale_colour_manual(values = grid.col,
                        name = '') +
    scale_size_continuous(range = c(2, 6),
                          limits = c(2, 40),
                          name = '-log10(p value)') +
    labs(x = 'r<sub>rm</sub>',
         y = '',
         subtitle = feature)
  
  return(p)
  
}
