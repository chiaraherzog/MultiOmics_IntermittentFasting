paired_change_compliance <- function(dat, variable = 'variable', nn = 4, p = 'p.format',
                                     colours = cols[c(4, 5, 1)],
                                     ylab = '',
                                     fdr = TRUE){
  complete <- dat |> 
    dplyr::filter(rowname == variable) |> 
    dplyr::group_by(subjectId) |> 
    dplyr::count() |> 
    ungroup() |> 
    dplyr::filter(n == nn)
  
  aspect.ratio = ifelse(nn == 4, 2, 1.5)
  
  if(fdr == T){
    # Load Wilcoxon test results with within-ome FDR
    load("out/wilcoxon-tests.Rdata")
    
    # Function
    get_fdr_compliance_pvals <- function(wilcoxon_tests = wilcoxon_tests){
        complist <- c('high', 'medium', 'low')
      
      out <- do.call(rbind, lapply(complist, function(comp){                # loop over compliance list
        
        long_pval <- wilcoxon_tests[[paste0(comp, ' compliance only')]] |>  # grab relevant wilcoxon test result item
          dplyr::filter(rowname == variable) |>                             # filter variable
          dplyr::select(rowname, ends_with("_adj")) |>                      # filter rowname and adjusted pvalue
          tidyr::pivot_longer(any_of(ends_with("_adj"))) |>                 # pivot to long format
          dplyr::rowwise() |> 
          dplyr::mutate(visitId = gsub("p_|_adj", "", name),                 # reformat p vlaues for plotting
                        compliance = paste0(comp, "<br>compliance"),
                        complabel = factor(compliance, levels = c("low<br>compliance", 'medium<br>compliance', 'high<br>compliance')),
                        p_value_label = ifelse(p == 'p.signif', gtools::stars.pval(value), signif(value, 3)),
                        p_value = value,
                        group1 = 'M0',
                        group2 = visitId) |> 
          dplyr::select(rowname, visitId, p_value, p_value_label, complabel, group1, group2) |> 
          dplyr::ungroup()
      }))
    
    }
  } 
  
  if(fdr == TRUE){
    
    # If FDR for high compliance only:
      fdr_p <- get_fdr_compliance_pvals(wilcoxon_tests, filter_high) |> 
        dplyr::filter(!is.na(p_value)) |> 
        dplyr::mutate(xmin = 1,
                      xmax = 2)
      
      # y position
      fdr_p$ypos <- max(dat_filtered$value + 0.05*dat_filtered$value)*0.95
      tip_lengths <- c(abs((unique(fdr_p$ypos) - max(dat_filtered$value)))*0.02,
                       abs((unique(fdr_p$ypos )- max(dat_filtered$value)))*0.02
      )
    
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
        ggpubr::stat_pvalue_manual(data = fdr_p,
                                   label = 'p_value_label',
                                   size = 2.7,
                                   tip.length = tip_lengths,
                                   bracket.size = 0.2, remove.bracket = F,
                                   y.position = 'ypos') +
        labs(x = "", y= ylab) +
        scale_colour_manual(values = colours,
                            name = "compliance",
                            aesthetics = c("fill", 'colour'))
    
  } else {
    
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
  } 
  
  return(p)
}
}