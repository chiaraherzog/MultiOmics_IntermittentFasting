plot_tissue_comparison <- function(dat_cerv,
                                   dat_buccal,
                                   dat_blood,
                                   variable = 'variable',
                                   variable_blood = NULL,
                                   nn = 4,
                                   colours = cols[c(6, 2, 1)],
                                   p = 'p.signif',
                                   ylab = '',
                                   filter_high = T,
                                   fdr = T){

  
  # Merge
  dat <- plyr::rbind.fill(dat_blood, dat_buccal, dat_cerv)  |> 
    dplyr::mutate(sampletype = gsub(".*: ", "", paste0(assay, "\nsample")))
  
  if(!is.null(variable_blood)){
    dat <- dat |> 
      dplyr::filter((sampletype %in% c('cervical\nsample', 'buccal\nsample') & rowname == variable) |
                      (sampletype == 'blood\nsample' & rowname == variable_blood))
  } else {
    dat <- dat |> 
      dplyr::filter(rowname == variable)
    }
  
  # Complete cases in each tissue
  complete_cerv <- dat_cerv |> 
    dplyr::filter(rowname == variable) |> 
    dplyr::group_by(subjectId) |> 
    dplyr::count() |> 
    ungroup() |> 
    dplyr::filter(n == nn)
  
  complete_buccal <- dat_buccal |> 
    dplyr::filter(rowname == variable) |> 
    dplyr::group_by(subjectId) |> 
    dplyr::count() |> 
    ungroup() |> 
    dplyr::filter(n == nn)
  
  complete_blood <- dat_blood |> 
    dplyr::filter(rowname == ifelse(!is.null(variable_blood), variable_blood, variable)) |> 
    dplyr::group_by(subjectId) |> 
    dplyr::count() |> 
    ungroup() |> 
    dplyr::filter(n == nn)
  
  stripcol <- strip_themed(background_x = elem_list_rect(fill = colours))
  
  if(filter_high == T){
    dat <- dat |>
      dplyr::filter(compliance == 'high')
  }
  
  dat <- dat |> dplyr::filter((sampletype == 'blood\nsample' & subjectId %in% complete_blood$subjectId) |
                                (sampletype == 'buccal\nsample' & subjectId %in% complete_buccal$subjectId) |
                                (sampletype == 'cervical\nsample' & subjectId %in% complete_cerv$subjectId)) |>
    dplyr::mutate(sampletype = factor(sampletype, levels = c('cervical\nsample', 'buccal\nsample', 'blood\nsample')))
  
  if(fdr){
    
    # Set up the variable filter
    filter <- if(!is.null(variable_blood)){
      "(assay %in% c('Composite methylation scores: buccal', 'Composite methylation scores: cervical') & x %in% variable) | (assay %in% c('Composite methylation scores: blood') & x %in% variable_blood)"
    } else {
      "(assay %in% c('Composite methylation scores: buccal', 'Composite methylation scores: cervical', 'Composite methylation scores: buccal') & x %in% variable)"
    }
    
    get_fdr_compliance_pvals <- function(wilcoxon_tests = wilcoxon_tests, filter_high = filter_high){
      
      if(filter_high == F){
        comp <- "high compliance only"
      } else {
        comp <- "overall"
      }
      
        long_pval <- wilcoxon_tests[[comp]] |>  # grab relevant wilcoxon test result item
          dplyr::filter(!!rlang::parse_expr(filter)) |>                               # filter variable
          dplyr::select(assay, x, ends_with("_adj")) |>                      # filter rowname and adjusted pvalue
          tidyr::pivot_longer(any_of(ends_with("_adj"))) |>                 # pivot to long format
          dplyr::rowwise() |> 
          dplyr::mutate(visitId = gsub("p_|_adj", "", name),  
                        p_value_label = ifelse(p == 'p.signif', gtools::stars.pval(value), signif(value, 3)),
                        p_value = paste0(signif(value, 2), gsub("[.]", "", p_value_label)),
                        group1 = 'M0',
                        group2 = visitId) |> 
          dplyr::select(x, assay, visitId, p_value, p_value_label, group1, group2) |> 
          dplyr::ungroup()
        
        return(long_pval)
    }
    
    fdr_p <- get_fdr_compliance_pvals(wilcoxon_tests, filter_high) |> 
      dplyr::filter(!is.na(p_value)) |> 
      dplyr::mutate(xmin = 1,
                    xmax = case_when(visitId == "M2" ~ 2,
                                     visitId == "M4" ~ 3,
                                     visitId == "M6" ~ 4))
    
    # y position
    ypos <- dat |>
      dplyr::group_by(assay) |> 
      dplyr::mutate(ypos = max(value),
                    maxval = max(value*0.02)) |>
      dplyr::ungroup() |>
      dplyr::select(assay, rowname, ypos, maxval, sampletype) |> dplyr::distinct()
    
    fdr_p <- fdr_p |> dplyr::left_join(ypos, by = c('x' = 'rowname',
                                                    'assay' = 'assay'))
    
    tip_lengths <- 0.02
    
    plot <- dat |>
      ggplot(aes(x = visitId,
                 y = value)) +
      geom_boxplot(outlier.shape = NA,
                   aes(fill = sampletype),
                   alpha = 0.2) +
      # geom_line(aes(group = subjectId,
      #               colour = sampletype),
      #           alpha = 0.2) +
      facet_wrap2(~sampletype,strip = stripcol) +
      scale_colour_manual(values = colours,
                          aesthetics = c('fill', 'colour')) +
      theme_bw() +
      theme(strip.text = element_text(face = 'bold',
                                      colour = 'white'),
            legend.position = 'none') +
      labs(x = '',
           y = ylab) +
      ggpubr::stat_pvalue_manual(
             data = fdr_p,
             hide.ns = T,
             bracket.nudge.y = c(rep(0, 3),
                                 rep(-0.45, 3),
                                 rep(-0.7, 3)),
             step.increase = 0.065,
             label = 'p_value',  # Column containing the p-value significance label
             xmin = 'xmin',            # x-axis start for the comparison bracket
             xmax = 'xmax',            # x-axis end for the comparison bracket
             y.position = 'ypos',      # y-position for the bracket
             size = 2.7,               # Font size for the p-value labels
             tip.length = 0.01,        # Adjust as necessary to control the length of the bracket tips
             bracket.size = 0.3,       # Adjust thickness of the brackets
             remove.bracket = FALSE    # Set to FALSE to keep the brackets visible  
             
           )
      
      
    } else {
  
  plot <- dat |> 
      ggplot(aes(x = visitId,
                 y = value)) +
      geom_boxplot(outlier.shape = NA,
                   aes(fill = sampletype),
                   alpha = 0.2) +
      geom_line(aes(group = subjectId,
                    colour = sampletype),
                alpha = 0.2) +
      facet_wrap2(~sampletype,
                  strip = stripcol) +
    scale_colour_manual(values = colours,
                        aesthetics = c('fill', 'colour')) +
      theme_bw() +
      theme(strip.text = element_text(face = 'bold',
                                      colour = 'white'),
            legend.position = 'none') +
      labs(x = '',
           y = ylab) +
      ggpubr::stat_compare_means(ref.group = 'M0',
                                 label = p,
                                 paired = T,
                                 hide.ns = T,
                                 label.y.npc = 0.95)
    } 
  
  return(plot)
    
}
