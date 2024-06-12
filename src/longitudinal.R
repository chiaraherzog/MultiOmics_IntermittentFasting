longitudinal <- function(dat,
                         variable = 'variable',
                         nn = 4,
                         p = 'p.format',
                         colour = 'grey',
                         ylab = ''){
  
  
  complete <- dat |> 
    dplyr::filter(rowname == variable) |> 
    dplyr::group_by(subjectId) |> 
    dplyr::count() |> 
    ungroup() |> 
    dplyr::filter(n == nn)
  
  
  if(nn == 2){
      plot <- dat |> 
        dplyr::filter(!is.na(compliance) & subjectId %in% complete$subjectId & rowname == variable) |> 
        ggplot(aes(x = visitId,
                   y = value)) +
        geom_boxplot(alpha = 0.3,
                     fill = colour,
                     outlier.shape = NA) +
        geom_line(aes(group = subjectId),
                  alpha = 0.03) +
        ggpubr::stat_compare_means(comparisons = list(c("M0", "M6")),
                                   label = p,
                                   size = 2.7,
                                   paired = T,
                                   label.y.npc = 0.95) +
        theme_bw() +
        theme(legend.position = 'none',
              axis.title.y = ggtext::element_markdown(),
              strip.text.x = ggtext::element_markdown(),
              aspect.ratio = 2) +
        labs(x = "", y = ylab) +
        scale_colour_manual(values = colour,
                            aesthetics = c("fill", 'colour'))
    } else {
      plot <- dat |> 
        dplyr::filter(!is.na(compliance) & subjectId %in% complete$subjectId & rowname == variable) |> 
        ggplot(aes(x = visitId,
                   y = value)) +
        geom_boxplot(alpha = 0.3,
                     fill = colour,
                     outlier.shape = NA) +
        geom_line(aes(group = subjectId),
                  alpha = 0.08) +
        ggpubr::stat_compare_means(ref.group = "M0",
                                   label = p,
                                   size = 2.7,
                                   paired = T,
                                   label.y.npc = 0.95) +
        theme_bw() +
        theme(legend.position = 'none',
              axis.title.y = ggtext::element_markdown(),
              strip.text.x = ggtext::element_markdown(),
              aspect.ratio = 2) +
        labs(x = "", y = ylab) +
        scale_colour_manual(values = colour,
                            aesthetics = c("fill", 'colour'))
    }
  
  return(plot)
}
