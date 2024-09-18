longitudinal <- function(dat,
                         variable = 'variable',
                         nn = 4,
                         filter_high = F,
                         p = 'p.format',
                         colour = 'grey',
                         ylab = '',
                         fdr = T){
  
  
  complete <- dat |> 
    dplyr::filter(rowname == variable) |> 
    dplyr::group_by(subjectId) |> 
    dplyr::count() |> 
    ungroup() |> 
    dplyr::filter(n == nn)
  
  if(fdr == T){
    # Load Wilcoxon test results with within-ome FDR
    load("out/wilcoxon-tests.Rdata")
    
    # Function
    get_fdr_compliance_pvals <- function(wilcoxon_tests = wilcoxon_tests, filter_high = filter_high){
      
      if(filter_high == F){
        complist <- 'overall'
        index = 'overall'
      } else {
        complist <- 'high compliance only'
      }
      
      out <- do.call(rbind, lapply(complist, function(comp){                # loop over compliance list
        
        long_pval <- wilcoxon_tests[[comp]] |>  # grab relevant wilcoxon test result item
          dplyr::filter(rowname == variable & assay %in% dat$assay) |>                             # filter variable
          dplyr::select(rowname, ends_with("_adj")) |>                      # filter rowname and adjusted pvalue
          tidyr::pivot_longer(any_of(ends_with("_adj"))) |>                 # pivot to long format
          dplyr::rowwise() |> 
          dplyr::mutate(visitId = gsub("p_|_adj", "", name),                 # reformat p vlaues for plotting
                        p_value_label = ifelse(p == 'p.signif', gtools::stars.pval(value),
                                               paste0(signif(value, 2), gsub("[.]", "", gtools::stars.pval(value)))),
                        p_value = paste0(signif(value), gsub("[.]", "", p_value_label)),
                        group1 = 'M0',
                        group2 = visitId) |> 
          dplyr::select(rowname, visitId, p_value, p_value_label, group1, group2) |> 
          dplyr::ungroup()
      }))
      
    }
    
  }
  
  
  if(nn == 2){
    
    dat_filtered <- dat |> 
      dplyr::filter(!is.na(compliance) & subjectId %in% complete$subjectId & rowname == variable)
      
    if(fdr == T){
      fdr_p <- get_fdr_compliance_pvals(wilcoxon_tests, filter_high) |> 
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
        geom_boxplot(alpha = 0.3,
                     fill = colour,
                     outlier.shape = NA) +
        geom_line(aes(group = subjectId),
                  alpha = 0.03) +
      ggpubr::stat_pvalue_manual(data = fdr_p,
                                 label = 'p_value_label',
                                 size = 2.6,
                                 tip.length = tip_lengths,
                                 bracket.size = 0.2, remove.bracket = F,
                                 y.position = 'ypos') +
        theme_bw() +
        theme(legend.position = 'none',
              axis.title.y = ggtext::element_markdown(),
              strip.text.x = ggtext::element_markdown(),
              aspect.ratio = 2) +
        labs(x = "", y = ylab) +
        scale_colour_manual(values = colour,
                            aesthetics = c("fill", 'colour'))
    } else {
      plot <- dat_filtered |> 
        ggplot(aes(x = visitId,
                   y = value)) +
        geom_boxplot(alpha = 0.3,
                     fill = colour,
                     outlier.shape = NA) +
        geom_line(aes(group = subjectId),
                  alpha = 0.03) +
        ggpubr::stat_compare_means(comparisons = list(c("M0", "M6")),
                                   label = p,
                                   size = 2.6,
                                   paired = T,
                                   label.y.npc = 0.95) +
        theme_bw() +
        theme(legend.position = 'none',
              axis.title.y = ggtext::element_markdown(),
              strip.text.x = ggtext::element_markdown(),
              aspect.ratio = 2) +
        labs(x = "", y = ylab) +
        scale_colour_manual(values = colour,
                            aesthetics = c("fill", 'colour'))
    }
  
    } else {
      
      dat_filtered <- dat |> 
        dplyr::filter(!is.na(compliance) & subjectId %in% complete$subjectId & rowname == variable) 
      
      if(fdr == T){
        fdr_p <- get_fdr_compliance_pvals(wilcoxon_tests, filter_high)
        
        # y position
        fdr_p$ypos <- max(dat_filtered$value + 0.05*dat_filtered$value)*0.95
      
      plot <- dat_filtered |> 
        ggplot(aes(x = visitId,
                   y = value)) +
        geom_boxplot(alpha = 0.3,
                     fill = colour,
                     outlier.shape = NA) +
        geom_line(aes(group = subjectId),
                  alpha = 0.08) +
        ggpubr::stat_pvalue_manual(data = fdr_p,
                                   label = 'p_value_label',
                                   x = 'visitId',
                                   
                                   size = 2.7,
                                   y.position = 'ypos') +
        theme_bw() +
        theme(legend.position = 'none',
              axis.title.y = ggtext::element_markdown(),
              strip.text.x = ggtext::element_markdown(),
              aspect.ratio = 2) +
        labs(x = "", y = ylab) +
        scale_colour_manual(values = colour,
                            aesthetics = c("fill", 'colour'))
    } else {
      plot <- dat_filtered |> 
        ggplot(aes(x = visitId,
                   y = value)) +
        geom_boxplot(alpha = 0.3,
                     fill = colour,
                     outlier.shape = NA) +
        geom_line(aes(group = subjectId),
                  alpha = 0.08) +
        ggpubr::stat_compare_means(ref.group = "M0",
                                   label = p,
                                   size = 2.7,
                                   paired = T,
                                   label.y.npc = 0.95) +
        theme_bw() +
        theme(legend.position = 'none',
              axis.title.y = ggtext::element_markdown(),
              strip.text.x = ggtext::element_markdown(),
              aspect.ratio = 2) +
        labs(x = "", y = ylab) +
        scale_colour_manual(values = colour,
                            aesthetics = c("fill", 'colour'))
    }
    }
  
  return(plot)
}
