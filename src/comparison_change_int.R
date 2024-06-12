comparison_change_int <- function(dat, variable = 'variable',
                                  ylab = '',
                                  p = "p.signif",
                                  colours = cols[c(4, 5, 1)],
                                  split_intervention = F){
    
    
    plot <- dat |> 
      dplyr::filter(visitId != "M0") |> 
      dplyr::filter(rowname == variable) |> 
      ggplot(aes(x = compliance,
                 y = value)) +
      geom_boxplot(alpha = 0.3,
                   aes(fill = compliance)) +
      ggbeeswarm::geom_beeswarm(aes(colour= compliance),
                                size = 0.9, alpha = 0.6) +
      ggpubr::stat_compare_means(comparisons = list(c("low", "high"),
                                                    c("low", "medium"),
                                                    c("medium", "high")),
                                 label = p,
                                 size = 2.7,
                                 method = 't.test',
                                 paired = F,
                                 hide.ns = T,
                                 label.y.npc = 0.95) +
      theme_bw() +
      theme(legend.position = 'none',
            axis.title.y = ggtext::element_markdown(),
            strip.text.x = ggtext::element_markdown(),
            aspect.ratio = 2) +
      labs(x = "", y = ylab) +
      scale_colour_manual(values = colours,
                          aesthetics = c("fill", 'colour')) +
      facet_wrap(~visitId)
    
    
    if(split_intervention==T){
      plot <- plot +
        facet_wrap(interventionId~visitId)
  
    } 
    
    return(plot)
    
  }
