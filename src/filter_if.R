filter_if <- function(data,
                      age_norm = F,
                      type_norm = F){
  
  d <- data |>
    dplyr::filter(interventionId!='S' & !visitId %in% c("M12", "M18")) |>
    dplyr::mutate(compliance = factor(compliance, levels = c("low", 'medium', 'high')))
  
  if(age_norm == T){
    d <- d |> 
    dplyr::mutate(age = case_when(visitId == "M0" ~ age_at_consent,
                                  visitId == "M2" ~ age_at_consent + 2/12,
                                  visitId == "M4" ~ age_at_consent + 4/12,
                                  visitId == "M6" ~ age_at_consent + 6/12))
    
    } 
  
  if(type_norm == T){
    d <- d |> 
      dplyr::mutate(type = ifelse(visitId == "M0", "Control",
                                "Follow-up"))
  }
  
  return(d)
}