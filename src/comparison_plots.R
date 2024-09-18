# Paired compliance
comparison_plots <- function(dat, variable = 'variable', nn = 4,
                                           p = 'p.format',
                                           colours = cols[c(2, 7)],
                                           ylab = '',
                             fdr = T,
                             sub = NULL){
  
  complete <- dat |>
    dplyr::filter(rowname == variable) |>
    dplyr::group_by(subjectId) |>
    dplyr::count() |>
    ungroup() |>
    dplyr::filter(n == nn)
  
  if(!is.null(sub)){
    dat <- dat |> 
      dplyr::filter(visitId == sub)
  }
  
  dat_filtered <- dat |> 
    dplyr::filter(grepl("high", compliance) & visitId != "M0" & subjectId %in% complete$subjectId) |> 
    dplyr::filter(rowname == variable)
  
  if(fdr == T){
    # Load Wilcoxon test results with within-ome FDR
    load("out/wilcoxon-tests.Rdata")
    
    # consider only clr transformed data microbiome ASVs, TSS transformed data microbiome families
    wilcoxon_tests <- lapply(wilcoxon_tests, function(df) {
      df %>% dplyr::filter(!assay %in% c("Saliva microbiome: families: clr","Saliva microbiome: ASVs",
                                         "Stool microbiome: families: clr","Stool microbiome: ASVs"))
    })
    
    get_fdr_compliance_pvals_IvK <- function(wilcoxon_tests = wilcoxon_tests){
      
      wilcoxon_tests[["I versus K"]] |>  # grab relevant wilcoxon test result item
          dplyr::filter(x == variable & grepl(unique(dat_filtered$assay), assay)) |>                             # filter variable
          dplyr::select(x, ends_with("_adj")) |>                      # filter rowname and adjusted pvalue
          tidyr::pivot_longer(any_of(ends_with("_adj"))) |>                 # pivot to long format
          dplyr::rowwise() |> 
          dplyr::mutate(visitId = gsub("p_|_adj", "", name),   
                        p_value_label = ifelse(p == 'p.signif', gtools::stars.pval(value), signif(value, 3)),
                        p_value = value,
                        group1 = 'I',
                        group2 = 'K') |> 
          dplyr::select(x, visitId, p_value, p_value_label, group1, group2) |> 
          dplyr::ungroup()
    }
    
    fdr_p <- get_fdr_compliance_pvals_IvK(wilcoxon_tests) |> 
      dplyr::filter(!is.na(p_value)) |> 
      dplyr::mutate(xmin = 1,
                    xmax = 2)
    
    # y position and brackets
    fdr_p$ypos <- max(dat_filtered$value + 0.05*dat_filtered$value)*0.95
    tip_lengths <- c(abs((unique(fdr_p$ypos) - max(dat_filtered$value)))*0.02,
                     abs((unique(fdr_p$ypos) - max(dat_filtered$value)))*0.02)
    
      
  plot <- dat |> 
    dplyr::filter(grepl("high", compliance) & visitId != "M0" & subjectId %in% complete$subjectId) |> 
    dplyr::filter(rowname == variable) |> 
    ggplot(aes(x = interventionId,
               y = value)) +
    geom_boxplot(alpha = 0.3,
                 aes(fill = interventionId),
                 outlier.shape = NA) +
    ggbeeswarm::geom_beeswarm(aes(colour= interventionId),
                              size = 0.9, alpha = 0.6) +
    ggpubr::stat_pvalue_manual(data = fdr_p,
                               label = 'p_value_label',
                               size = 2.7,
                               tip.length = tip_lengths,
                               bracket.size = 0.2, remove.bracket = F,
                               y.position = 'ypos') +
    theme_bw() +
    theme(legend.position = 'none',
          axis.title.y = ggtext::element_markdown(),
          strip.text.x = ggtext::element_markdown(),
          aspect.ratio = 2
          ) +
    labs(x = "", y = ylab) +
    scale_colour_manual(values = colours,
                        aesthetics = c("fill", 'colour'))+
    facet_wrap(~visitId)
  
  
  } else {
    
  
  plot <- dat |> 
      dplyr::filter(grepl("high", compliance) & visitId != "M0" & subjectId %in% complete$subjectId) |> 
      dplyr::filter(rowname == variable) |> 
      ggplot(aes(x = interventionId,
                 y = value)) +
      geom_boxplot(alpha = 0.3,
                   aes(fill = interventionId),
                   outlier.shape = NA) +
      ggbeeswarm::geom_beeswarm(aes(colour= interventionId),
                                size = 0.9, alpha = 0.6) +
      ggpubr::stat_compare_means(comparisons = list(c("I", "K")),
                                 label = p,
                                 size = 2.7,
                                 method = 'wilcox',
                                 paired = F,
                                 label.y.npc = 0.95) +
      theme_bw() +
      theme(legend.position = 'none',
            axis.title.y = ggtext::element_markdown(),
            strip.text.x = ggtext::element_markdown(),
            aspect.ratio = 2) +
      labs(x = "", y = ylab) +
      scale_colour_manual(values = colours,
                          aesthetics = c("fill", 'colour'))+
      facet_wrap(~visitId)
  
  }
  
  return(plot)
  
}
