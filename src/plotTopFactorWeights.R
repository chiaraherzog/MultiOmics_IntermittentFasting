#' plotTopFactorWeights
#' 
#' @name plotTopFactorWeights
#' @param mofa mofa object
#' @param factorFilter  which factor to filter - numeric
#' @param orientation portrait or landscape?
#' @param top_n how many top features to include
#' @param col_low colour for low values (negative)
#' @param col_mid colour for 0
#' @param col_high colour for high values (positive)
#'
#' @returns a list object with the plot and data (top weights)
plotTopFactorWeights <- function(mofa, factorFilter,
                                 orientation = 'portrait',
                                 top_n = 30,
                                 col_low = '#71b5a9',
                                 col_mid = 'white',
                                 col_high = '#bd647d'){
  
  # Extract weights
  w <- get_weights(mofa,  as.data.frame = T, scale = T) |>
    dplyr::mutate(absval = abs(value)) |> 
    dplyr::filter(factor %in% paste0("Factor", factorFilter))
  
  # Filter & fix labels
  load("src/vars.Rdata")
  
  top_w <- w |> 
    dplyr::arrange(absval) |> 
    dplyr::slice_max(absval, n = top_n) |> 
    dplyr::ungroup() |> 
    dplyr::mutate(feature = gsub("_Composite methylation scores: buccal|_Composite methylation scores: cervical|_Composite methylation scores: blood|_Saliva microbiome: families|_Urine nuclear magnetic resonance: normalized|_Saliva nuclear magnetic resonance: normalized|_Methylation PCs: cervical|_Methylation PCs: blood|_Methylation PCs: buccal|_Stool microbiome: families", "", feature)) |> 
    dplyr::left_join(vars, by = c("feature" = "x",
                                  "view" = "assay")) |> 
    dplyr::mutate(aname2 = case_when(grepl("methylation", view) & grepl("buccal", view) ~ "Buccal methylation scores",
                                     grepl("methylation", view) & grepl("blood", view) ~ "Blood methylation scores",
                                     grepl("methylation", view) & grepl("cervical", view) ~ "Cervical methylation scores",
                                     grepl("Saliva", view, ignore.case = T) & grepl("nuclear", view) ~ "Saliva metabolome",
                                     grepl("Saliva", view, ignore.case = T) & !grepl("metabolome", view) ~ "Saliva microbiome",
                                     grepl("Urine", view, ignore.case = T) & grepl("nuclear", view) ~ "Urine metabolome",
                                     grepl("Stool", view, ignore.case = T) ~ "Faecal microbiome",
                                     TRUE ~ view),
                  label = ifelse(is.na(label), feature, label),
                  label = ifelse(grepl("Flow", aname2), gsub("[+]", "<sup>+</sup>", label),
                                 label),
                  label = ifelse(grepl("Flow", aname2), gsub("[-]", "<sup>-</sup>", label),
                                 label)
    )
  
  
  if(orientation == 'landscape'){
    
    plot <- top_w |> 
      dplyr::arrange(absval) |> 
      ggplot(aes(x = forcats::fct_reorder(label, -absval),
                 y = forcats::fct_reorder(aname2, absval),
                 fill = value)) +
      geom_tile(color = 'black') +
      theme_bw() +
      theme(axis.text.x = element_markdown(angle = 90,
                                           hjust = 1,
                                           vjust = 0.5,
                                           color = 'black'),
            axis.text.y = element_markdown(hjust = 1),
            axis.ticks.y = element_blank(),
            axis.ticks.x = element_blank(),
            legend.position = 'top',
            legend.key.height = unit(4, "mm"),
            panel.grid = element_blank()) +
      scale_fill_gradient2(low = col_low,
                           mid = col_mid,
                           high = col_high,
                           midpoint = 0,
                           limits = c(-1, 1)
      ) +
      coord_cartesian(expand = F) +
      labs(x = '',
           y = '')
    
  } else {
    
    plot <- top_w |> 
      dplyr::arrange(absval) |> 
      ggplot(aes(y = forcats::fct_reorder(label, absval),
                 x = forcats::fct_reorder(aname2, absval),
                 fill = value)) +
      geom_tile(color = 'black') +
      theme_bw() +
      theme(axis.text.x = element_markdown(angle = 90,
                                           hjust = 1,
                                           vjust = 0.5,
                                           color = 'black'),
            axis.text.y = element_markdown(hjust = 1),
            axis.ticks.y = element_blank(),
            axis.ticks.x = element_blank(),
            legend.position = 'top',
            legend.key.height = unit(4, "mm"),
            panel.grid = element_blank()) +
      scale_y_discrete(position = 'right') +
      scale_fill_gradient2(low = col_low,
                           mid = col_mid,
                           high = col_high,
                           midpoint = 0,
                           limits = c(-1, 1)
      ) +
      coord_cartesian(expand = F) +
      labs(x = '',
           y = '')
    
  }
  
  plot <- plot + coord_fixed()
  
  return(list(plot = plot,
              data = top_w))
}
