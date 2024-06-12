lmm <- function(dat, variable,
                colours = cols[c(4, 5, 1)],
                ylab = '',
                nn = 4){
  
  library(lmerTest)
  library(lme4)
  
  complete <- dat |> 
    dplyr::filter(rowname == variable) |> 
    dplyr::group_by(subjectId) |> 
    dplyr::count() |> 
    ungroup() |> 
    dplyr::filter(n == nn)
  
  mod <- lmer(value ~ time*compliance + age_at_consent + (1 | subjectId), data = dat[dat$rowname==variable & dat$subjectId %in% complete$subjectId,])
  
  # Plot
  if(nn == 4){
    
  p <- ggeffects::ggpredict(mod, terms = c("time [2, 4, 6]", "compliance [low, medium, high]"))
  
  plot <- plot(p) +
    theme_bw() +
    theme(legend.position = 'none',
          axis.title.y = ggtext::element_markdown()) +
    labs(x = "", y = ylab,
         title = "") +
    scale_colour_manual(values = colours,
                        aesthetics = c("fill", 'colour'),
                        name = 'Compliance group') +
    scale_x_continuous(breaks = c(2, 4, 6),
                       labels = c("M2", "M4", "M6"))
  
  } else {
    p <- ggeffects::ggpredict(mod, terms = c("time [0,6]", "compliance [low, medium, high]"))
    
    plot <- plot(p) +
      theme_bw() +
      theme(legend.position = 'none',
            axis.title.y = ggtext::element_markdown()) +
      labs(x = "", y = ylab,
           title = "") +
      scale_colour_manual(values = colours,
                          aesthetics = c("fill", 'colour'),
                          name = 'Compliance group') +
      scale_x_continuous(breaks = c(6),
                         labels = c("M6"))
  }
  
  # Out
  return(list(plot = plot,
              model = mod))
}