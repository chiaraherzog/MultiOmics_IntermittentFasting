#' @name comparison_plots_menopause
#' @description
#' A function to plot pre versus postmenopausal women (high compliance)
#' @param dat data frame in longFormat
#' @param variable rowname variable to use
#' @param nn Number of observations (if only measure twice, it is 2, if 4 measurements, 4). 4 by default
#' @param p Input for ggpubr p label - p.format (default) or p.signif (stars)
#' @param colours Colours for pre and postmenopausal
#' @param ylab Y label axis text (can be richtext for ggtext)
#' @param subsetM optional - should any specific visit be subset?
#' @param all_months Should M0 be included or not? (Default is true but may not be needed if showing only change scores.)
#' 
comparison_plots_menopause <- function(dat, variable = 'variable', nn = 4,
                             p = 'p.format',
                             colours = c("#5f70a8", "#bd647d"),
                             ylab = '',
                             subsetM = NULL,
                             all_months = T){
  
  
  if(!"mpstatrs" %in% colnames(dat)){
    stop("column mpstatrs not present. please check the input data")
  }
  
  
  complete <- dat |>
    dplyr::filter(rowname == variable) |>
    dplyr::group_by(subjectId) |>
    dplyr::count() |>
    ungroup() |>
    dplyr::filter(n == nn)
  
  if(!is.null(subsetM)){
    dat <- dat |> 
      dplyr::filter(visitId == subsetM)
  }
  
  if(all_months == F){
    dat <- dat |> 
    dplyr::filter(grepl("high", compliance) & visitId != "M0" & subjectId %in% complete$subjectId)
  } else {
    dat <- dat |> 
      dplyr::filter(grepl("high", compliance) & subjectId %in% complete$subjectId)
  }
  
  plot <- dat |> 
    dplyr::mutate(mpstatrs = ifelse(mpstatrs == 'yes', 'post', 'pre'),
                  mpstatrs = factor(mpstatrs, levels = c('pre', 'post'))) |> 
    dplyr::filter(rowname == variable) |> 
    ggplot(aes(x = mpstatrs,
               y = value)) +
    geom_boxplot(alpha = 0.3,
                 aes(fill = mpstatrs),
                 outlier.shape = NA) +
    ggbeeswarm::geom_beeswarm(aes(colour= mpstatrs),
                              size = 0.9, alpha = 0.6) +
    ggpubr::stat_compare_means(comparisons = list(c("pre", "post")),
                               label = p,
                               size = 2.7,
                               method = 'wilcox',
                               paired = F,
                               label.y.npc = 0.85) +
    theme_bw() +
    theme(legend.position = 'none',
          axis.title.y = ggtext::element_markdown(),
          strip.text.x = ggtext::element_markdown(),
          aspect.ratio = 2) +
    labs(x = "", y = ylab) +
    scale_colour_manual(values = colours,
                        aesthetics = c("fill", 'colour'))+
    facet_wrap(~visitId, nrow = 1) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))
  
  return(plot)
  
}
