
plotFP <- function(data, variable, nn = 4,
                   p = 'p.format',
                   colours = c('gold2', 'gold4', 'deepskyblue4', 'deepskyblue3', 'grey'),
                   ylab = '',
                   filter_high = F){
  
  if(filter_high == T){
    data <- dat |> dplyr::filter(compliance == 'high')
  }
  
  complete <- data |> 
    dplyr::filter(rowname == variable) |> 
    dplyr::group_by(subjectId) |> 
    dplyr::count() |> 
    ungroup() |> 
    dplyr::filter(n == nn)
  
  p <-  data |> dplyr::filter(subjectId %in% complete$subjectId & rowname == variable) |> 
    dplyr::mutate(pattern = ifelse(is.na(pattern), "unknown", pattern),
                  pattern = stringr::str_wrap(pattern, width = 10),
                  pattern = factor(pattern, levels = c("breakfast\ncancellation",
                                                       "alternating\n(predominant\nbreakfast)",
                                                       "dinner\ncancellation",
                                                       "alternating\n(predominant\ndinner)",
                                                       "unknown"))
                  ) |> 
    ggplot(aes(x = visitId,
               y = value)) +
    geom_boxplot(aes(fill = pattern),
                 alpha = 0.2) +
    geom_line(aes(group = subjectId,
                  colour = pattern)) +
    facet_grid(interventionId ~ pattern) +
    ggpubr::stat_compare_means(paired = T, label = 'p.format',
                               label.y.npc = 0.95,
                               size = 2.4) +
    theme_bw() +
    theme_bw() +
    theme(legend.position = 'none',
          axis.title.y = ggtext::element_markdown(),
          strip.text = element_text(size = 6.4)
          )  +
    labs(x = "", y = ylab) +
    scale_colour_manual(values = colours,
                        aesthetics = c('colour', 'fill'))
  
  return(p)
                               
  }