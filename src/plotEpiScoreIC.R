#' @param df dataframe with methylation values to be used for plotting. can be df_scaled or df_raw
#' @param variable epigenetic variable to be plotted on the y axis
#' @param byCompliance parameter indicating whether compliance groups should be split. T by default
#' @param complete parameter indicating whether only complete observations should be used
#' @param ylim optional feature indicating y lims (coord cartesian)
#' @param colours to be used for visitIds
#' @param ylab optional feature indicating ylab. using ggtext/element markdown so can include html formatting.
#' @returns plot with epigenetic score (y) versus ic (x)


plotEpiScoreIC <- function(df,
                           variable,
                           byCompliance = T,
                           complete = T,
                           ylim = NULL,
                           colours = c("#48a0af",
                                       "#1b69a1",
                                       "#ec6669",
                                       "#f39668"),
                           iVk = F,
                           ylab = NULL){
  
  tmp <- df |> 
    dplyr::filter(rowname %in% c(variable, "ic")) |> 
    tidyr::pivot_wider(id_cols = c("subjectId", "visitId", "compliance", "interventionId"),
                       names_from = "rowname",
                       values_from = "value")
  
  # complete
  if(complete == T){
    comp <- tmp |> 
      dplyr::group_by(subjectId) |> 
      dplyr::filter(all(c("M0", "M2", "M4", "M6") %in% visitId)) |> 
      dplyr::pull(subjectId) |> 
      unique()
    
    tmp <- tmp |> 
      dplyr::filter(subjectId %in% comp)
  
    }
  
  if(is.null(ylab)){ylab <- variable}
  
  if(iVk == T){
    
    p <- tmp |> 
      dplyr::filter(compliance == 'high') |> 
      ggplot(aes(x = ic,
                 y = .data[[variable]],
                 colour = visitId)) +
      geom_point(alpha = 0.5) +
      geom_smooth(method = 'lm',
                  se = F) +
      facet_wrap(~interventionId) +
      coord_cartesian(ylim = ylim) +
      theme_bw() +
      theme(axis.text.y = element_markdown()) +
      coord_cartesian(ylim = ylim) +
      scale_colour_manual(values = colours,
                          name = 'visit') +
      labs(x = 'immune cell proportion',
           ylab = ylab)
    
  } else if(byCompliance == T){
    
   p <- tmp |> 
      ggplot(aes(x = ic,
                 y = .data[[variable]],
                 colour = visitId)) +
      geom_point(alpha = 0.5) +
      geom_smooth(method = 'lm',
                  se = F) +
      facet_wrap(~compliance) +
      coord_cartesian(ylim = ylim) +
      theme_bw() +
      theme(axis.text.y = element_markdown()) +
      coord_cartesian(ylim = ylim) +
      scale_colour_manual(values = colours,
                          name = 'visit') +
      labs(x = 'immune cell proportion',
           ylab = ylab)
    
  } else {
    
   p <- tmp |> 
      ggplot(aes(x = ic,
                 y = .data[[variable]],
                 colour = visitId)) +
      geom_point(alpha = 0.5) +
      geom_smooth(method = 'lm',
                  se = F) +
      theme_bw() +
      theme(axis.text.y = element_markdown()) +
      coord_cartesian(ylim = ylim) +
      scale_colour_manual(values = colours,
                          name = 'visit') +
      labs(x = 'immune cell proportion',
           ylab = ylab)
  
    }
  
  return(p)
  
}
