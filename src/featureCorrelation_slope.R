featureCorrelation_slope <- function(slope,
                               x, y,
                               corMethod = 'pearson'){
  
  tmp <- slope |> 
    tibble::rownames_to_column('feature')   |> 
    dplyr::filter(grepl(paste0(c(x, y), collapse = "|"), feature, ignore.case = T)) |> 
    dplyr::mutate(feature = ifelse(grepl(x, feature, ignore.case = T), x,
                                   y)) |> 
    dplyr::mutate(feature = gsub(",| |-|\\$|^", "", feature)) |> 
    tibble::column_to_rownames('feature') |> 
    t() |> 
    as.data.frame()
  
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
