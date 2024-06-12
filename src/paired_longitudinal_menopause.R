#' @name paired_longitudinal_menopause
#' @description
#'A function to plot changes in pre and post menopausal women separately
#' @param dat data frame in longFormat
#' @param variable rowname variable to use
#' @param nn Number of observations (if only measure twice, it is 2, if 4 measurements, 4). 4 by default
#' @param p Input for ggpubr p label - p.format (default) or p.signif (stars)
#' @param colours Colours for pre and postmenopausal
#' @param ylab Y label axis text (can be richtext for ggtext)
#' @param filter_high Should high compliance be filtered? T by default.
paired_longitudinal_menopause <- function(dat, variable = 'variable', nn = 4,
                                           p = 'p.format',
                                           colours = c("#5f70a8", "#bd647d"),
                                           ylab = '',
                                           filter_high = F){
  
  if(!"mpstatrs" %in% colnames(dat)){
    stop("column mpstatrs not present. please check the input data")
  }
  
  complete <- dat |> 
    dplyr::filter(rowname == variable) |> 
    dplyr::filter(!is.na(value)) |> 
    dplyr::group_by(subjectId) |> 
    dplyr::count() |> 
    ungroup() |> 
    dplyr::filter(n == nn)
  
  
  if(filter_high == T){
    dat <- dat |> 
      dplyr::filter(grepl("high", compliance) & subjectId %in% complete$subjectId) |> 
      dplyr::filter(rowname == variable)
  } else {
    dat <- dat |> 
      dplyr::filter(subjectId %in% complete$subjectId)|> 
      dplyr::filter(rowname == variable)
  }
  
  plot <- dat |> 
      dplyr::mutate(mpstatrs = ifelse(mpstatrs == 'yes', 'post', 'pre'),
                    mpstatrs = factor(mpstatrs, levels = c('pre', 'post'))) |> 
      dplyr::filter(rowname == variable) |>
      ggplot(aes(x = visitId,
                 y = value)) +
      geom_line(aes(group = subjectId,
                    colour = mpstatrs),
                alpha = 0.3) +
      geom_boxplot(alpha = 0.3,
                   aes(fill = mpstatrs),
                   outlier.shape = NA) +
      ggpubr::stat_compare_means(ref.group = "M0",
                                 label = "p.signif",
                                 size = 2.7,
                                 paired = T,
                                 label.y.npc = 0.95) +
    labs(x = "", y = ylab) +  
    theme_bw() +
      theme(legend.position = 'none',
            aspect.ratio = 2,
            axis.title.y = ggtext::element_markdown(),
            strip.text.x = ggtext::element_markdown()) +
      scale_colour_manual(values = colours,
                          aesthetics = c("fill", 'colour')) +
      facet_wrap(~mpstatrs) +
      scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))
  
  
  return(plot)
  
}
