
plotFP <- function(data, variable, nn = 4,
                   p = 'p.format',
                   colours = c('gold2', 'gold4', 'deepskyblue4', 'deepskyblue3', 'grey'),
                   ylab = '',
                   filter_high = F,
                   fdr = TRUE){
  
  if(filter_high == T){
    data <- dat |> dplyr::filter(compliance == 'high')
  }
  
  complete <- data |> 
    dplyr::filter(rowname == variable) |> 
    dplyr::group_by(subjectId) |> 
    dplyr::count() |> 
    ungroup() |> 
    dplyr::filter(n == nn)
  
  if(fdr == T){
    # Load Wilcoxon test results with within-ome FDR
    load("out/wilcoxon-tests.Rdata")
    
    # Function
    get_fdr_compliance_pvals <- function(wilcoxon_tests = wilcoxon_tests){
      
      list <- c("Breakfast cancellation",
                "Dinner cancellation",
                "Alternating (predominant breakfast cancellation)",
                "Alternating (predominant dinner cancellation)",
                "Fasting pattern unknown")
      
      out <- do.call(rbind, lapply(list, function(comp){                # loop over compliance list
        
        long_pval <- wilcoxon_tests[[comp]] |>  # grab relevant wilcoxon test result item
          dplyr::filter(rowname == variable) |>                             # filter variable
          dplyr::select(rowname, ends_with("_adj")) |>                      # filter rowname and adjusted pvalue
          tidyr::pivot_longer(any_of(ends_with("_adj"))) |>                 # pivot to long format
          dplyr::rowwise() |> 
          dplyr::mutate(visitId = gsub("p_|_adj", "", name),                 # reformat p vlaues for plotting
                        pattern = gsub(" cancellation)", ")", comp),
                        pattern = gsub("Fasting pattern ", "", pattern),
                        pattern = tolower(pattern),
                        pattern = stringr::str_wrap(pattern, width = 10),
                        pattern = factor(pattern, levels = c("breakfast\ncancellation",
                                                             "alternating\n(predominant\nbreakfast)",
                                                             "dinner\ncancellation",
                                                             "alternating\n(predominant\ndinner)",
                                                             "unknown")),
                        p_value_label = ifelse(p == 'p.signif', gtools::stars.pval(value), signif(value, 3)),
                        p_value = value,
                        group1 = 'M0',
                        group2 = visitId) |> 
          dplyr::select(rowname, visitId, p_value, p_value_label, pattern, group1, group2) |> 
          dplyr::ungroup()
      }))
      
    }
    
  }
  
  dat_filtered <- data |>
    dplyr::filter(subjectId %in% complete$subjectId & rowname == variable) |> 
    dplyr::mutate(pattern = ifelse(is.na(pattern), "unknown", pattern),
                  pattern = stringr::str_wrap(pattern, width = 10),
                  pattern = factor(pattern, levels = c("breakfast\ncancellation",
                                                       "alternating\n(predominant\nbreakfast)",
                                                       "dinner\ncancellation",
                                                       "alternating\n(predominant\ndinner)",
                                                       "unknown")))
  
  if(fdr == T){
    
    fdr_p <- get_fdr_compliance_pvals(wilcoxon_tests) |> 
      dplyr::filter(!is.na(p_value)) |> 
      dplyr::mutate(xmin = 1,
                    xmax = 2)
    
    # y position
    fdr_p$ypos <- max(dat_filtered$value + 0.05*dat_filtered$value)*0.95
    tip_lengths <- c(abs((unique(fdr_p$ypos) - max(dat_filtered$value)))*0.02,
                     abs((unique(fdr_p$ypos )- max(dat_filtered$value)))*0.02
    )
  
  plot <- dat_filtered |> 
    ggplot(aes(x = visitId,
               y = value)) +
    geom_boxplot(aes(fill = pattern),
                 alpha = 0.2) +
    geom_line(aes(group = subjectId,
                  colour = pattern)) +
    facet_grid(~ pattern) +
    ggpubr::stat_pvalue_manual(data = fdr_p,
                               label = 'p_value_label',
                               size = 2.7,
                               tip.length = tip_lengths,
                               bracket.size = 0.2, remove.bracket = F,
                               y.position = 'ypos') +
    theme_bw() +
    theme_bw() +
    theme(legend.position = 'none',
          axis.title.y = ggtext::element_markdown(),
          strip.text = element_text(size = 6.4)
          )  +
    labs(x = "", y = ylab)  +
  
    scale_colour_manual(values = colours,
                        aesthetics = c('colour', 'fill'))
  
                     
  } else {
    
    plot <- dat_filtered |> 
      ggplot(aes(x = visitId,
                 y = value)) +
      geom_boxplot(aes(fill = pattern),
                   alpha = 0.2) +
      geom_line(aes(group = subjectId,
                    colour = pattern)) +
      facet_grid(~ pattern) +
      ggpubr::stat_compare_means(comparisons = list(c("M0", "M6")),
                                 size = 2.7, 
                                 paired = T,
                                 label = p) +
      theme_bw() +
      theme_bw() +
      theme(legend.position = 'none',
            axis.title.y = ggtext::element_markdown(),
            strip.text = element_text(size = 6.4)
      )  +
      labs(x = "", y = ylab)  +
      scale_colour_manual(values = colours,
                          aesthetics = c('colour', 'fill'))
  }
  
  return(plot)
}
  
