# Paired compliance
paired_intervention_compliance <- function(dat, variable = 'variable', nn = 4,
                                           p = 'p.format',
                                           colours = cols[c(2, 8)],
                                           ylab = ''){
  
  complete <- dat |> 
    dplyr::filter(rowname == variable) |> 
    dplyr::group_by(subjectId) |> 
    dplyr::count() |> 
    ungroup() |> 
    dplyr::filter(n == nn)
  
  if(nn == 4){
  plot <- dat |> 
      dplyr::filter(compliance == 'high' & subjectId %in% complete$subjectId & rowname == variable) |>
      dplyr::mutate(complabel = paste0(compliance, "<br>compliance"),
                    complabel = factor(complabel, levels = c("low<br>compliance", 'medium<br>compliance', 'high<br>compliance'))) |> 
      ggplot(aes(x = visitId,
                 y = value)) +
      geom_boxplot(alpha = 0.3,
                   aes(fill = interventionId)) +
      geom_line(aes(group = subjectId,
                    colour = interventionId),
                alpha = 0.3) +
      ggpubr::stat_compare_means(ref.group = "M0",
                                 label = p,
                                 size = 2.7,
                                 paired = T,
                                 label.y.npc = 0.95) +
      theme_bw() +
      theme(legend.position = 'none',
            aspect.ratio = 2,
            axis.title.y = ggtext::element_markdown(),
            strip.text.x = ggtext::element_markdown()) +
      labs(x = "", y = ylab) +
      scale_colour_manual(values = colours,
                          aesthetics = c("fill", 'colour')) +
      facet_wrap(~interventionId)
    } else if(nn == 2){
      plot <- dat |> 
        dplyr::filter(grepl("high", compliance)) |> 
        dplyr::filter(!is.na(compliance) & subjectId %in% complete$subjectId & rowname == variable) |> 
        ggplot(aes(x = visitId,
                   y = value)) +
        geom_boxplot(alpha = 0.3,
                     aes(fill = compliance)) +
        geom_line(aes(group = subjectId,
                      colour = compliance),
                  alpha = 0.3) +
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
        scale_colour_manual(values = colours[3],
                            aesthetics = c("fill", 'colour'))
  }
  
  return(plot)
  
}
