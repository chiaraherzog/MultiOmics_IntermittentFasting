featureCorrelation_absolute <- function(data,
                                        experiments,
                                        visit,
                                        x, y,
                                        corMethod = 'pearson'){
  
  tmp <- longFormat(data[,data$visitId == visit & data$interventionId != 'S', experiments]) |> 
    as.data.frame() |> 
    dplyr::mutate(feature = paste0(assay, "_", rowname)) |> 
    dplyr::filter(grepl(paste0(c(x, y), collapse = "|"), feature, ignore.case = T)) |> 
    dplyr::mutate(feature = ifelse(grepl(x, feature, ignore.case = T), x,
                                   y)) |> 
    dplyr::mutate(feature = gsub(",| |-|\\$|^", "", feature)) |> 
    tidyr::pivot_wider(id_cols = 'primary',
                       names_from = 'feature',
                       values_from = 'value')
  
  x = gsub(",| |-|\\$|^", "", x)
  y <- gsub(",| |-|\\$|^", "", y)
  
  p <- tmp |> 
    ggplot(aes_string(x = x,
                      y = y)) +
    geom_point() +
    geom_smooth(method = 'lm', se = F) +
    ggpubr::stat_cor(method = corMethod)
  
  print(p)
}
