# Paired compliance
paired_longitudinal_compliance <- function(dat, variable = 'variable', nn = 4,
                                           p = 'p.format',
                                           colours = cols[c(4, 5, 1)],
                                           ylab = '',
                                           filter_high = F,
                                           fdr = TRUE){
  
  complete <- dat |> 
    dplyr::filter(rowname == variable) |> 
    dplyr::group_by(subjectId) |> 
    dplyr::count() |> 
    ungroup() |> 
    dplyr::filter(n == nn)
  
  if(fdr == T){
    # Load Wilcoxon test results with within-ome FDR
    load("out/wilcoxon-tests.Rdata")
    
    # consider only clr transformed data microbiome ASVs, TSS transformed data microbiome families
    wilcoxon_tests <- lapply(wilcoxon_tests, function(df) {
      df %>% dplyr::filter(!assay %in% c("Saliva microbiome: families: clr","Saliva microbiome: ASVs",
                                         "Stool microbiome: families: clr","Stool microbiome: ASVs"))
    })
    
    # Function
    get_fdr_compliance_pvals <- function(wilcoxon_tests = wilcoxon_tests, filter_high = filter_high){
      
      if(filter_high == F){
        complist <- c('high', 'medium', 'low')
      } else {
        complist <- 'high'
      }
      
      out <- do.call(rbind, lapply(complist, function(comp){                # loop over compliance list
        
        long_pval <- wilcoxon_tests[[paste0(comp, ' compliance only')]] |>  # grab relevant wilcoxon test result item
          #dplyr::filter(rowname == variable) |>                             # filter variable
          dplyr::filter(x == variable) |> 
          #dplyr::select(rowname, ends_with("_adj")) |>                      # filter rowname and adjusted pvalue
          dplyr::select(x, ends_with("_adj")) |> 
          tidyr::pivot_longer(any_of(ends_with("_adj"))) |>                 # pivot to long format
          dplyr::rowwise() |>
          dplyr::mutate(visitId = gsub("p_|_adj", "", name),                 # reformat p vlaues for plotting
                        compliance = paste0(comp, "<br>compliance"),
                        complabel = factor(compliance, levels = c("low<br>compliance", 'medium<br>compliance', 'high<br>compliance')),
                        p_value_label = ifelse(p == 'p.signif', gtools::stars.pval(value), signif(value, 3)),
                        p_value = value,
                        group1 = 'M0',
                        group2 = visitId) |>
          #dplyr::select(rowname, visitId, p_value, p_value_label, complabel, group1, group2) |>
          dplyr::select(x, visitId, p_value, p_value_label, complabel, group1, group2) %>%
          dplyr::ungroup()
      }))
      
    }
  
    }
    
  # If filter high is false, keep all complete cases and facet by compliance
  
  if(filter_high == F){
    
    if(nn == 2){
      
      dat_filtered <- dat |> 
        dplyr::filter(!is.na(compliance) & subjectId %in% complete$subjectId & rowname == variable) |>
        dplyr::mutate(complabel = paste0(compliance, "<br>compliance"),
                      complabel = factor(complabel, levels = c("low<br>compliance", 'medium<br>compliance', 'high<br>compliance')))
      
      # If FDR for high compliance only:
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
          geom_line(aes(group = subjectId,
                        colour = compliance),
                    alpha = 0.3) +
          geom_boxplot(alpha = 0.3,
                       aes(fill = compliance),
                       outlier.shape = NA) +
          theme_bw() +
          theme(legend.position = 'none',
                aspect.ratio = 2,
                axis.title.y = ggtext::element_markdown(),
                strip.text.x = ggtext::element_markdown()) +
          labs(x = "", y = ylab) +
          scale_colour_manual(values = colours,
                              aesthetics = c("fill", 'colour')) +
          ggpubr::stat_pvalue_manual(data = fdr_p,
                                     label = 'p_value_label',
                                     size = 2.7,
                                     tip.length = tip_lengths,
                                     bracket.size = 0.2, remove.bracket = F,
                                     y.position = 'ypos') +
          facet_wrap(~complabel)
        
      } else {
        
        plot <- dat_filtered |> 
          ggplot(aes(x = visitId,
                     y = value)) +
          geom_line(aes(group = subjectId,
                        colour = compliance),
                    alpha = 0.3) +
          geom_boxplot(alpha = 0.3,
                       aes(fill = compliance),
                       outlier.shape = NA) +
          theme_bw() +
          theme(legend.position = 'none',
                aspect.ratio = 2,
                axis.title.y = ggtext::element_markdown(),
                strip.text.x = ggtext::element_markdown()) +
          labs(x = "", y = ylab) +
          scale_colour_manual(values = colours,
                              aesthetics = c("fill", 'colour')) +
          ggpubr::stat_compare_means(ref.group = "M0",
                                     label = p,
                                     size = 2.7,
                                     paired = T,
                                     label.y.npc = 0.95) +
          facet_wrap(~complabel)
      }
      
      
    } else {
    
    # filter compliance and set up labels
    dat_filtered <- dat |>
      dplyr::filter(!is.na(compliance) & subjectId %in% complete$subjectId & rowname == variable) |>
      dplyr::mutate(complabel = paste0(compliance, "<br>compliance"),
                    complabel = factor(complabel, levels = c("low<br>compliance", 'medium<br>compliance', 'high<br>compliance')))
    
    # Get FDR p values if fdr = T ----
    # p values: FDR corrected or not
    if(fdr == T){
      fdr_p <- get_fdr_compliance_pvals(wilcoxon_tests, filter_high)
      # y position
      fdr_p$ypos <- max(dat_filtered$value + 0.04*dat_filtered$value)*0.95
      
      plot <- dat_filtered |> 
        ggplot(aes(x = visitId,
                   y = value)) +
        geom_line(aes(group = subjectId,
                      colour = compliance),
                  alpha = 0.3) +
        geom_boxplot(alpha = 0.3,
                     aes(fill = compliance),
                     outlier.shape = NA) +
        theme_bw() +
        theme(legend.position = 'none',
              aspect.ratio = 2,
              axis.title.y = ggtext::element_markdown(),
              strip.text.x = ggtext::element_markdown()) +
        labs(x = "", y = ylab) +
        scale_colour_manual(values = colours,
                            aesthetics = c("fill", 'colour')) +
        ggpubr::stat_pvalue_manual(data = fdr_p,
                                   label = 'p_value_label',
                                   x = 'visitId',
                                   size = 2.7,
                                   y.position = 'ypos') +
        facet_wrap(~complabel)
      
      } else {
      
      plot <- dat_filtered |> 
        ggplot(aes(x = visitId,
                   y = value)) +
        geom_line(aes(group = subjectId,
                      colour = compliance),
                  alpha = 0.3) +
        geom_boxplot(alpha = 0.3,
                     aes(fill = compliance),
                     outlier.shape = NA) +
        theme_bw() +
        theme(legend.position = 'none',
              aspect.ratio = 2,
              axis.title.y = ggtext::element_markdown(),
              strip.text.x = ggtext::element_markdown()) +
        labs(x = "", y = ylab) +
        scale_colour_manual(values = colours,
                            aesthetics = c("fill", 'colour')) +
        ggpubr::stat_compare_means(ref.group = "M0",
                                   label = p,
                                   size = 2.7,
                                   paired = T,
                                   label.y.npc = 0.95) +
        facet_wrap(~complabel)
      
    }
   
    }
  } else {
    
    if(nn == 2){
      
      dat_filtered <- dat |> 
        dplyr::filter(!is.na(compliance) & subjectId %in% complete$subjectId & rowname == variable & grepl("high", compliance)) 
      
      # If FDR for high compliance only:
      if(fdr == T){
        fdr_p <- get_fdr_compliance_pvals(wilcoxon_tests, filter_high) |> 
          dplyr::filter(!is.na(p_value)) |> 
          dplyr::mutate(xmin = 1,
                        xmax = 2)
        
        # y position
        fdr_p$ypos <- max(dat_filtered$value + 0.06*dat_filtered$value)*0.95
        tip_lengths <- c(abs((fdr_p$ypos - max(dat_filtered[dat_filtered$visitId=='M0',]$value)))*0.01,
                         abs((fdr_p$ypos - max(dat_filtered[dat_filtered$visitId=='M6',]$value)))*0.01)
        
        plot <- dat_filtered |> 
          ggplot(aes(x = visitId,
                     y = value)) +
          geom_line(aes(group = subjectId,
                        colour = compliance),
                    alpha = 0.3) +
          geom_boxplot(alpha = 0.3,
                       aes(fill = compliance)) +
          ggpubr::stat_pvalue_manual(data = fdr_p,
                                     label = 'p_value_label',
                                     size = 2.7,
                                     bracket.size = 0.2, remove.bracket = F,
                                     tip.length = tip_lengths,
                                     y.position = 'ypos') +
          theme_bw() +
          theme(legend.position = 'none',
                axis.title.y = ggtext::element_markdown(),
                strip.text.x = ggtext::element_markdown(),
                aspect.ratio = 2) +
          labs(x = "", y = ylab) +
          scale_colour_manual(values = colours[3],
                              aesthetics = c("fill", 'colour'))
        
      } else {
          
      plot <- dat |> 
        dplyr::filter(!is.na(compliance) & subjectId %in% complete$subjectId & rowname == variable) |> 
        ggplot(aes(x = visitId,
                   y = value)) +
        geom_line(aes(group = subjectId,
                      colour = compliance),
                  alpha = 0.3) +
        geom_boxplot(alpha = 0.3,
                     aes(fill = compliance)) +
        ggpubr::stat_compare_means(comparisons = list(c("M0", "M6")),
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
        scale_colour_manual(values = colours[3],
                            aesthetics = c("fill", 'colour'))
      }
      
    } else {
      
      dat_filtered <- dat |> 
        dplyr::filter(grepl("high", compliance)) |> 
        dplyr::filter(!is.na(compliance) & subjectId %in% complete$subjectId & rowname == variable)
      
      # If FDR for high compliance only:
      if(fdr == T){
        fdr_p <- get_fdr_compliance_pvals(wilcoxon_tests, filter_high)
        # y position
        fdr_p$ypos <- max(dat_filtered$value + 0.04*dat_filtered$value)*0.95
        
        plot <- dat_filtered |>
          ggplot(aes(x = visitId,
                     y = value)) +
          geom_line(aes(group = subjectId,
                        colour = compliance),
                    alpha = 0.3) +
          geom_boxplot(alpha = 0.3,
                       aes(fill = compliance),
                       outlier.shape = NA) +
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
          scale_colour_manual(values = colours[3],
                              aesthetics = c("fill", 'colour'))
        
        } else {
    
        plot <- dat_filtered |> 
          ggplot(aes(x = visitId,
                     y = value)) +
          geom_line(aes(group = subjectId,
                        colour = compliance),
                    alpha = 0.3) +
          geom_boxplot(alpha = 0.3,
                       aes(fill = compliance),
                       outlier.shape = NA) +
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
          scale_colour_manual(values = colours[3],
                              aesthetics = c("fill", 'colour'))
        
      }
    
   }
  }
  return(plot)
}
