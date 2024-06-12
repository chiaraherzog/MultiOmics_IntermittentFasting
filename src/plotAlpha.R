## Helper function for plotting alpha-diversity over time
# Author: Charlotte Vavourakis

plotAlpha <- function(dat, rarmin, category){
  
  require(tidyverse)
  require(ampvis2)
  require(ggpubr)
  require(patchwork)
  
  alphadiversity = amp_alphadiv(dat, measure = c("shannon", "simpson","observed","invsimpson"), rarefy = rarmin, richness = TRUE)
  
  # make long format
  out = alphadiversity %>%
    dplyr::select(subjectId,visitId, Shannon, Chao1) %>%
    pivot_longer(
      cols = c("Shannon","Chao1"),
      names_to = "index",
      values_to = "score",
      values_drop_na = TRUE
    )
  
  compare = list( c("M0", "M2"), c("M0", "M4"), c("M0", "M6") )
  
  ## richness
  df = out %>%
    filter(index=="Chao1") %>%
    dplyr::select(-index) %>%
    droplevels() 
  
  df$visitId = factor(df$visitId, levels=c("M0","M2","M4","M6"))
  
  plot1 = ggboxplot(df, x = "visitId", y = "score") + 
    geom_line(aes(group = subjectId), colour = "#1b69a1", alpha = 0.3) + 
    stat_compare_means(comparisons = compare, paired=TRUE)+ 
    stat_compare_means(label.y = 50) +
    ylab("Chao1 (richness)") +
    xlab("") +
    theme(legend.position="none") +
    ggtitle(category)
  
  ## diversity
  df = out %>%
    filter(index=="Shannon") %>%
    dplyr::select(-index) %>%
    droplevels() 
  
  df$visitId = factor(df$visitId, levels=c("M0","M2","M4","M6"))
  
  plot2 = ggboxplot(df, x = "visitId", y = "score") + 
    geom_line(aes(group = subjectId), colour = "#1b69a1", alpha = 0.3) + 
    stat_compare_means(comparisons = compare, paired=TRUE)+ 
    stat_compare_means(label.y = 0) +
    ylab("Shannon (diversity)") +
    xlab("") +
    theme(legend.position="none") +
    ggtitle("\n")
  
  ## combine
  
  #plot = plot1 + plot2 + plot_annotation(title = category)
  plot = plot1 + plot2
  
  return(plot)
  
}