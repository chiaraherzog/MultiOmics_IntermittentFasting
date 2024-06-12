#' @param m1 Measure 1
#' @param m2 Measure 2
#' @param lab1 Label for measure 1
#' @param lab2 Label for measure 2
#' @return ggplot with repeated measures correlation

rmcorr_plot <- function(m1, m2, lab1, lab2){
  
  library(ggtext)
  # Data are loaded
  load(here("out/plot_rmcorr_examples.Rdata"))
  dat <- dat |> 
    dplyr::select(subjectId, visitId, any_of(m1), any_of(m2)) |> 
    dplyr::filter(!is.na(.data[[m1]]) & !is.na(.data[[m2]]))
  
  # Colours
  load(here("src/cols_for_assays.Rdata"))

  # run rmcorr
  mod <- rmcorr::rmcorr(participant = subjectId,
                        measure1 = m1, measure2 = m2,
                        dataset = dat)
  
# Rename labels
  labels <- data.frame(id = c(1, 2),
                       measures = c(m1, m2),
                       name = c(lab1, lab2))
  
  labels <- labels |> 
    tidyr::separate(measures, "_", into = c('assay', 'featureid'), remove = F, extra = 'merge') |> 
    dplyr::mutate(assay = case_when(
      grepl("haemogram", assay) | (grepl("exam", assay) & !grepl("_fe|_fv|_sysbp|_diabp|_vo2max|_rel|_abs|_max", featureid)) ~ "Routine bloods",
      grepl("cytometry", assay) ~ "Flow cytometry: immune cells",
      grepl("magnetic", assay) & grepl("Saliva", assay) ~ "Saliva metabolome",
      grepl("magnetic", assay) & grepl("Urine", assay) ~ "Urine metabolome",
      grepl("methylation", assay) & grepl("cervical", assay) ~ "Cervical methylation",
      grepl("methylation", assay) & grepl("buccal", assay) ~ "Buccal methylation",
      grepl("methylation", assay) & grepl("blood", assay) ~ "Blood methylation",
      grepl("microbiome", assay) & grepl("Saliva", assay) ~ "Saliva microbiome",
      grepl("microbiome", assay) & grepl("Stool", assay) ~ "Stool microbiome",
      grepl("composition", assay) ~ "Body composition",
      (grepl("exam", assay) & grepl("_fe|_fv|_sysbp|_diabp|_vo2max|_rel|_abs|_max", featureid)) | (grepl("sono", assay) & grepl("pwv|imt|plaque", featureid)) ~ "Functional clinical features")) 
  
  # Add colour
  labels <- labels |> 
    dplyr::mutate(assay = paste0("<span style='color:",
                                 cols_for_assays[match(assay, names(cols_for_assays))],
                                 "'>", assay, "</span>"))
  
  # Full label
  labels <- labels |> 
    dplyr::mutate(label = paste0("<b>", name, "</b><br>(", assay, ")"))
  
  r <- paste0("r<sub><i>rm</i></sub>", " = ",
              format(mod$r, digits = 2),
              " (p=", format(mod$p, digits = 2), ")")
  
  # plot
  plot <- dat |> 
    dplyr::filter(!is.na(.data[[m1]] & !is.na(.data[[m2]]))) |> 
    dplyr::select(subjectId, visitId, .data[[m1]], .data[[m2]]) |> 
    ggplot(aes(x = .data[[m1]],
               y = .data[[m2]])) +
    geom_point(aes(colour = subjectId),
               size = 1.5,
               alpha = 0.9) +
    geom_line(aes(y = mod$model$fitted.values,
                  colour = subjectId),
              linetype = 1,
              alpha = 0.5) +
    geom_smooth(method = 'lm', se = F,
                linetype = 'dotted',
                colour = 'black') +
    scale_colour_manual(values = grDevices::colorRampPalette(cols[c(8, 1,2,3,5,4,6,7)])(114),
                        aesthetics = c('colour', 'fill')) +
    theme_bw() +
    theme(legend.position = 'none',
          axis.title = element_markdown()) +
    labs(x = labels$label[1],
         y = labels$label[2]) +
    annotate('richtext', label = r,
             Inf, -Inf,
             fill = NA, label.color = NA,
             hjust = 2.2,
             vjust = -12)

  return(plot)
  }
