plot_tissue_comparison <- function(dat_cerv,
                                   dat_buccal,
                                   dat_blood,
                                   variable = 'variable',
                                   variable_blood = NULL,
                                   nn = 4,
                                   colours = cols[c(6, 2, 1)],
                                   p = 'p.signif',
                                   ylab = '',
                                   filter_high = T){
  
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
  
  plot <- dat |> dplyr::filter((sampletype == 'blood\nsample' & subjectId %in% complete_blood$subjectId) |
                           (sampletype == 'buccal\nsample' & subjectId %in% complete_buccal$subjectId) |
                           (sampletype == 'cervical\nsample' & subjectId %in% complete_cerv$subjectId)) |>
      dplyr::mutate(sampletype = factor(sampletype, levels = c('cervical\nsample', 'buccal\nsample', 'blood\nsample'))) |> 
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
  
  return(plot)
    
}
