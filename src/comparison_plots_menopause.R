comparison_plots_menopause <- function(dat, variable = 'variable', nn = 4,
                                       p = 'p.format',
                                       colours = c("#5f70a8", "#bd647d"),
                                       ylab = '',
                                       subsetM = NULL,
                                       all_months = T,
                                       fdr = FALSE){  # Added fdr argument
  
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
  
  dat <- dat |> 
    dplyr::mutate(mpstatrs = ifelse(mpstatrs == 'yes', 'post', 'pre'),
                  mpstatrs = factor(mpstatrs, levels = c('pre', 'post'))) |> 
    dplyr::filter(rowname == variable)
  
  if(fdr == TRUE){
    # Load Wilcoxon test results with FDR correction
    load("out/wilcoxon-tests.Rdata")
    
    get_fdr_pvals_PreM_PostM <- function(wilcoxon_tests = wilcoxon_tests){
      
      wilcoxon_tests[["PreM versus PostM"]] |>  # grab relevant Wilcoxon test result item
        dplyr::filter(x == variable & grepl(unique(dat$assay), assay)) |>  # filter for the variable
        dplyr::select(x, ends_with("_adj")) |>  # select adjusted p-value
        tidyr::pivot_longer(any_of(ends_with("_adj"))) |>  # pivot to long format
        dplyr::rowwise() |> 
        dplyr::mutate(visitId = gsub("p_|_adj|_yes", "", name),
                      p_value_label = ifelse(p == 'p.signif', gtools::stars.pval(value), signif(value, 3)),
                      p_value = value,
                      group1 = 'PreM',
                      group2 = 'PostM') |> 
        dplyr::select(x, visitId, p_value, p_value_label, group1, group2) |> 
        dplyr::ungroup()
    }
    
    fdr_p <- get_fdr_pvals_PreM_PostM(wilcoxon_tests) |> 
      dplyr::filter(!is.na(p_value)) |> 
      dplyr::mutate(xmin = 1, xmax = 2)
    
    # y position and brackets
    fdr_p$ypos <- max(dat$value + 0.05 * dat$value) * 0.95
    tip_lengths <- c(abs((unique(fdr_p$ypos) - max(dat$value)))*0.02,
                     abs((unique(fdr_p$ypos) - max(dat$value)))*0.02)
    
    plot <- ggplot(dat, aes(x = mpstatrs, y = value)) +
      geom_boxplot(alpha = 0.3, aes(fill = mpstatrs), outlier.shape = NA) +
      ggbeeswarm::geom_beeswarm(aes(colour = mpstatrs), size = 0.9, alpha = 0.6) +
      ggpubr::stat_pvalue_manual(data = fdr_p,
                                 label = 'p_value_label',
                                 size = 2.7,
                                 tip.length = tip_lengths,
                                 bracket.size = 0.2,
                                 remove.bracket = F,
                                 y.position = 'ypos') +
      theme_bw() +
      theme(legend.position = 'none',
            axis.title.y = ggtext::element_markdown(),
            strip.text.x = ggtext::element_markdown(),
            aspect.ratio = 2) +
      labs(x = "", y = ylab) +
      scale_colour_manual(values = colours,
                          aesthetics = c("fill", 'colour')) +
      facet_wrap(~visitId, nrow = 1) +
      scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))
    
  } else {
    plot <- ggplot(dat, aes(x = mpstatrs, y = value)) +
      geom_boxplot(alpha = 0.3, aes(fill = mpstatrs), outlier.shape = NA) +
      ggbeeswarm::geom_beeswarm(aes(colour = mpstatrs), size = 0.9, alpha = 0.6) +
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
                          aesthetics = c("fill", 'colour')) +
      facet_wrap(~visitId, nrow = 1) +
      scale_y_continuous(expand = expansion(mult = c(0.05, 0.15)))
  }
  
  return(plot)
}
