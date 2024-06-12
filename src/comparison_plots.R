# Paired compliance
comparison_plots <- function(dat, variable = 'variable', nn = 4,
                                           p = 'p.format',
                                           colours = cols[c(2, 7)],
                                           ylab = '',
                             sub = NULL){
  
  complete <- dat |>
    dplyr::filter(rowname == variable) |>
    dplyr::group_by(subjectId) |>
    dplyr::count() |>
    ungroup() |>
    dplyr::filter(n == nn)
  
  if(!is.null(sub)){
    dat <- dat |> 
      dplyr::filter(visitId == sub)
  }
  
    plot <- dat |> 
      dplyr::filter(grepl("high", compliance) & visitId != "M0" & subjectId %in% complete$subjectId) |> 
      dplyr::filter(rowname == variable) |> 
      ggplot(aes(x = interventionId,
                 y = value)) +
      geom_boxplot(alpha = 0.3,
                   aes(fill = interventionId),
                   outlier.shape = NA) +
      ggbeeswarm::geom_beeswarm(aes(colour= interventionId),
                                size = 0.9, alpha = 0.6) +
      ggpubr::stat_compare_means(comparisons = list(c("I", "K")),
                                 label = p,
                                 size = 2.7,
                                 method = 'wilcox',
                                 paired = F,
                                 label.y.npc = 0.95) +
      theme_bw() +
      theme(legend.position = 'none',
            axis.title.y = ggtext::element_markdown(),
            strip.text.x = ggtext::element_markdown(),
            aspect.ratio = 2) +
      labs(x = "", y = ylab) +
      scale_colour_manual(values = colours,
                          aesthetics = c("fill", 'colour'))+
      facet_wrap(~visitId)
  
  return(plot)
  
}
