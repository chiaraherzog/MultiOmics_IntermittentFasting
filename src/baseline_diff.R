baseline_diff <- function(dat, variables = ''){
  
  # Compute differences
  
  for(variable in variables){
  
    if(variable %in% c("weight", "bmi", "fm", "bcm", "ecw",
                       'leukocytes', 'neutrophils', 'erythrocytes', 'hematocrit', 'hemoglobin', 'mch', 'mchc', 'mcv',
                       'distrwidth', 'platelets', 'monocytes', 'granim', 'crp')){
      nn = 4
    } else {
      nn = 2
    }
    
    complete <- dat |> 
      dplyr::group_by(subjectId) |> 
      dplyr::reframe(n = sum(!is.na(get(UQ(variable))))) |> 
      ungroup() |> 
      dplyr::filter(n == nn)
    
    tmp1 <- dat |> 
      tidyr::pivot_wider(id_cols = c(interventionId, subjectId, compliance_cat),
                         names_from = visitId,
                         values_from = {{variable}}) |> 
      dplyr::mutate(across(any_of(c("M2", "M4", "M6")), ~ . - M0))
    
    tab1_all <- tmp1 |> 
      dplyr::reframe(across(any_of(c("M0", "M2", "M4", "M6")),
                              ~ paste0(round(mean(., na.rm = T), 2),
                                       " (",
                                       round(sd(., na.rm = T), 2), ")")),
                       variable = paste0(variable)) |> 
      dplyr::relocate(variable)
    
    # Complete cases:
    tab2_complete <- tmp1 |> 
      dplyr::filter(subjectId %in% complete$subjectId & !is.na(compliance_cat)) |> 
      dplyr::reframe(across(any_of(c("M0", "M2", "M4", "M6")),
                            ~ paste0(round(mean(., na.rm = T), 2),
                                     " (",
                                     round(sd(., na.rm = T), 2), ")")),
                       variable = paste0(variable)) |> 
      dplyr::relocate(variable)
    
    # Split by compliance group
    tab3_compliance <- tmp1 |> 
      dplyr::filter(subjectId %in% complete$subjectId & !is.na(compliance_cat)) |> 
      dplyr::group_by(compliance_cat) |> 
      dplyr::reframe(across(any_of(c("M0", "M2", "M4", "M6")),
                            ~ paste0(round(mean(., na.rm = T), 2),
                                     " (",
                                     round(sd(., na.rm = T), 2), ")")),
                       variable = paste0(variable),
                       compliance_cat = compliance_cat) |> 
      distinct() |> 
      dplyr::relocate(variable, compliance_cat)
    
    # Split by intervention
    tab4_intervention <- tmp1 |> 
      dplyr::filter(subjectId %in% complete$subjectId & !is.na(compliance_cat)) |> 
      dplyr::group_by(interventionId) |> 
      dplyr::reframe(across(any_of(c("M0", "M2", "M4", "M6")),
                            ~ paste0(round(mean(., na.rm = T), 2),
                                     " (",
                                     round(sd(., na.rm = T), 2), ")")),
                     variable = paste0(variable),
                     interventionId = interventionId) |> 
      distinct() |> 
      dplyr::relocate(variable, interventionId)
    
    # Split by intervention and compliance
    tab5_detail <- tmp1 |> 
      dplyr::filter(subjectId %in% complete$subjectId & !is.na(compliance_cat)) |> 
      dplyr::group_by(interventionId, compliance_cat) |> 
      dplyr::reframe(across(any_of(c("M0", "M2", "M4", "M6")),
                            ~ paste0(round(mean(., na.rm = T), 2),
                                     " (",
                                     round(sd(., na.rm = T), 2), ")")),
                     variable = paste0(variable),
                     interventionId = interventionId,
                     compliance_cat = compliance_cat) |> 
      distinct() |> 
      dplyr::relocate(variable, interventionId, compliance_cat)
    
    if(variable == variables[1]){
      tab1 <- tab1_all
      tab2 <- tab2_complete
      tab3 <- tab3_compliance
      tab4 <- tab4_intervention
      tab5 <- tab5_detail
    } else {
      tab1 <- plyr::rbind.fill(tab1, tab1_all)
      tab2 <- plyr::rbind.fill(tab2, tab2_complete)
      tab3 <- plyr::rbind.fill(tab3, tab3_compliance)
      tab4 <- plyr::rbind.fill(tab4, tab4_intervention)
      tab5 <- plyr::rbind.fill(tab5, tab5_detail)
    }
    
    }
  
  return(list(tab1 = tab1,
              tab2 = tab2, 
              tab3 = tab3,
              tab4 = tab4,
              tab5 = tab5))
    
}
