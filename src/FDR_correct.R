
# helper function for p-value correction
# input is dataframe with lmm results in long format with columns estimate and p_value

FDR_correct <- function(pdat_long){
  
  # Separate the data so that each term is corrected separately (factors also separately)
  dfs <- list()
  
  for (i in 1:length(unique(pdat_long$estimate))){
    dfs[[i]] = pdat_long %>% filter(estimate == unique(pdat_long$estimate)[i])
    
  }
  
  # Apply FDR
  fdr_correction <- function(data) {
    p_values = data$p_value
    adjusted_p_values = p.adjust(p_values, method = "fdr")
    data$p_value = adjusted_p_values
    return(data)
  }
  
  dfs <- lapply(dfs,fdr_correction)
  
  # combine dfs
  pdat_long_corrected = do.call(rbind,dfs)
  
  return(pdat_long_corrected)
}

# alternative implementation, combining categories for each term
# FDR_correct <- function(pdat_long){
#   
#   # Separate the data into the groups that should be corrected separately
#   df_age = pdat_long %>% filter(grepl("age", estimate, ignore.case = TRUE))
#   df_bmi = pdat_long %>% filter(grepl("bmi", estimate, ignore.case = TRUE))
#   df_visitId = pdat_long %>% 
#     filter(!grepl("age|bmi|:", estimate, ignore.case = TRUE)) %>%
#     filter(grepl("visitId", estimate, ignore.case = TRUE))
#   df_extra_term = pdat_long %>% 
#     filter(!grepl("age|bmi|visitId|:", estimate, ignore.case = TRUE))
#   df_interaction = pdat_long %>% 
#     filter(!grepl("age|bmi", estimate, ignore.case = TRUE)) %>%
#     filter(grepl(":", estimate, ignore.case = TRUE))
#   
#   # Apply FDR
#   fdr_correction <- function(data) {
#     p_values = data$p_value
#     adjusted_p_values = p.adjust(p_values, method = "fdr")
#     data$p_value = adjusted_p_values
#     return(data)
#   }
#   
#   df_age = fdr_correction(df_age)
#   df_bmi = fdr_correction(df_bmi)
#   df_visitId = fdr_correction(df_visitId)
#   df_extra_term = fdr_correction(df_extra_term)
#   df_interaction = fdr_correction(df_interaction)
#   
#   # Combine the corrected data frames
#   pdat_long_corrected = bind_rows(df_age, df_bmi, df_visitId, df_extra_term, df_interaction)
#   
#   return(pdat_long_corrected)
# }