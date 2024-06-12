paired_change_compliance <- function(dat, variable = 'variable', nn = 4, p = 'p.format',
                                     colours = cols[c(4, 5, 1)],
                                     ylab = ''){
  complete <- dat |> 
    dplyr::filter(rowname == variable) |> 
    dplyr::group_by(subjectId) |> 
    dplyr::count() |> 
    ungroup() |> 
    dplyr::filter(n == nn)
  
  aspect.ratio = ifelse(nn == 4, 2, 1.5)
  
  p <- dat |>
    dplyr::filter(subjectId %in% complete$subjectId & rowname == variable & visitId != 'M0') |> 
    ggplot(aes(x = compliance,
               y = value,
               fill = compliance)) +
    geom_boxplot(alpha = 0.5) +
    ggbeeswarm::geom_quasirandom(aes(colour = compliance),
                                 size = 0.5,
                                 alpha = 0.7) +
    facet_wrap(~visitId) +
    geom_hline(yintercept = 0,
               linetype = 'dotted') +
    theme_bw() +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.y = element_markdown(),
          aspect.ratio = aspect.ratio) +
    ggpubr::stat_compare_means(size = 2.7,
                               comparisons = list(c("high", "medium"),
                                                  c("high", "low"),
                                                  c("low", "medium"))) +
    labs(x = "", y= ylab) +
    scale_colour_manual(values = colours,
                        name = "compliance",
                        aesthetics = c("fill", 'colour'))
  
  return(p)
}
