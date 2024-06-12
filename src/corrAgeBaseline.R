corrAgeBaseline <- function(data){
  
  # Keep only those with n>=3
  filter <- data |> 
    dplyr::filter(visitId == "M0") |> 
    dplyr::group_by(rowname) |> 
    dplyr::reframe(n = sum(!is.na(value)),
                   sd = sd(value, na.rm = T)) |> 
    dplyr::filter(n >= 3 & sd != 0) 
  
  cat("Continuing with ", nrow(filter), " variables (SD > 0 and n>=3).", sep = "")
  
  tmp <- data |> 
    dplyr::filter(rowname %in% filter$rowname & visitId == "M0") |> 
    dplyr::group_by(rowname) |> 
    dplyr::reframe(cor = cor(value, age_at_consent),
                   p = cor.test(value, age_at_consent)$p.value) |> 
    dplyr::ungroup()
  
  return(tmp)
  
}
