
#' @param exp which experiment should be plotted
#' @param lmm_data High compliance intervention output. (out_lmm$`Intervention (higher compliance)`)
#' @param cols Colour palette (leave blank for default colours [recommended], but can be altered if desired)
#' @param relabel (optional). Data frame supplying labels for the plot. x = variable, label = relabelled variable
#' @param age_cor should correlation with age be plotted on the side? default to F
#' @param bmi_cor should correlation with bmi be plotted on the side? default to F
#' @param FDR should p-values be corrected? default to F
#' @return ComplexHeatmap object
#' 
# version made for microbiomes
# option to plot age and/or bmi correlations at baseline
# option to FDR correct p-values, for each term (and category) separately 
# all terms with p.corrected < 0.05 are plotted.
# for uncorrected p-values, all rows with p < 0.01 are plotted

combination_heatmap_high_factor <- function(exp,
                                     lmm_data,
                                     cols = NULL,
                                     relabel = NULL,
                                     age_cor = F,
                                     bmi_cor = F,
                                     FDR = F){
  
  source("src/FDR_correct.R")
  
  if (is.null(cols)){
    cols <- c("#1b69a1", "#48a0af", "#71b5a9", "#ec6669", "#f39668", "#bd647d", "#832c9b", "#5f70a8")
  }
  
  ## load output lmm models and correct p-values if requested, keep only significant results ##----
  
  if (FDR == F) {
  
    tmp <- lmm_data |> 
      filter(str_detect(x, exp)) %>%
      dplyr::select(x, 
                    `estimate_visitIdM2`, `p.value_visitIdM2`,
                    `estimate_visitIdM2:interventionIdK`, `p.value_visitIdM2:interventionIdK`,
                    `estimate_visitIdM4`, `p.value_visitIdM4`,
                    `estimate_visitIdM4:interventionIdK`, `p.value_visitIdM4:interventionIdK`,
                    `estimate_visitIdM6`, `p.value_visitIdM6`,
                    `estimate_visitIdM6:interventionIdK`, `p.value_visitIdM6:interventionIdK`) 
    
    tmp$x = gsub(paste0(exp,"_"),"",tmp$x)
                    
    # filter out only significant rows and replace pval with * symbols
    tmp <- tmp %>%
      filter(rowSums(across(starts_with("p.value") & (contains("visitId")|contains("interventionId")), ~ . < 0.01)) > 0) %>% 
      dplyr::mutate(across(contains("p.value"), ~ ifelse(is.na(.), "", gtools::stars.pval(.)))) |> 
      dplyr::mutate(across(contains("p.value"), ~ gsub("[.]", "", .)))
    
  } else if (FDR == T){
    
    tmp <- lmm_data |> 
      filter(str_detect(x, exp))
    
    tmp_p <- tmp %>%
      dplyr::select(x,contains("p.value")) %>%
      pivot_longer(cols = contains("p.value"),names_to = "estimate", values_to = "p_value") %>%
      FDR_correct() %>%
      pivot_wider(id_cols = x, names_from = estimate, values_from = p_value)
    
    tmp <- tmp %>%
      select(!contains("p.value")) %>%
      left_join(tmp_p, by = "x") %>%dplyr::select(x, 
                                                  `estimate_visitIdM2`, `p.value_visitIdM2`,
                                                  `estimate_visitIdM2:interventionIdK`, `p.value_visitIdM2:interventionIdK`,
                                                  `estimate_visitIdM4`, `p.value_visitIdM4`,
                                                  `estimate_visitIdM4:interventionIdK`, `p.value_visitIdM4:interventionIdK`,
                                                  `estimate_visitIdM6`, `p.value_visitIdM6`,
                                                  `estimate_visitIdM6:interventionIdK`, `p.value_visitIdM6:interventionIdK`) 
    
    tmp$x = gsub(paste0(exp,"_"),"",tmp$x)
    
    # filter out only significant rows and replace pval with * symbols
    tmp <- tmp %>%
      filter(rowSums(across(starts_with("p.value") & (contains("visitId")|contains("interventionId")), ~ . < 0.05)) > 0) %>% 
      dplyr::mutate(across(contains("p.value"), ~ ifelse(is.na(.), "", gtools::stars.pval(.)))) |> 
      dplyr::mutate(across(contains("p.value"), ~ gsub("[.]", "", .)))
    
  }
  
  ## baseline correlations ##----
  
  if(age_cor == T & bmi_cor == F){
    load("out/corrAgeBaseline.R")
    corr <- corr |> 
      tidyr::separate(rowname, "_", into = c("assay", "name"), extra = 'merge') |> 
      dplyr::slice(match(tmp$x, name))
    
    tmp <- tmp |> 
      dplyr::left_join(dplyr::select(corr, name, cor, p), by = c('x' = 'name'))
    
    hr <- rowAnnotation('Correlation with age\n(baseline)' = anno_points(tmp$cor,
                                                                         pch = ifelse(tmp$p < 0.05, 19, 1),
                                                                         ylim = c(-1, 1),
                                                                         gp = gpar(col = ifelse(tmp$cor<0, cols[1], cols[4]))))
  } else if (age_cor == F & bmi_cor == T) {
    load("out/corrBmiBaseline.R")
    corr <- corr |> 
      tidyr::separate(rowname, "_", into = c("assay", "name"), extra = 'merge') |> 
      dplyr::slice(match(tmp$x, name))
    
    tmp <- tmp |> 
      dplyr::left_join(dplyr::select(corr, name, cor, p), by = c('x' = 'name'))
    
    hr <- rowAnnotation('Correlation with bmi\n(baseline)' = anno_points(tmp$cor,
                                                                         pch = ifelse(tmp$p < 0.05, 19, 1),
                                                                         ylim = c(-1, 1),
                                                                         gp = gpar(col = ifelse(tmp$cor<0, cols[1], cols[4]))))
  } else if (age_cor == T & bmi_cor == T) {
    load("out/corrAgeBaseline.R")
    corr <- corr |> 
      tidyr::separate(rowname, "_", into = c("assay", "name"), extra = 'merge') |> 
      dplyr::slice(match(tmp$x, name)) %>%
      dplyr::rename(corr_age = cor) %>%
      dplyr::rename(p_age = p)
    
    tmp <- tmp |> 
      dplyr::left_join(dplyr::select(corr, name, corr_age, p_age), by = c('x' = 'name'))
    
    load("out/corrBmiBaseline.R")
    corr <- corr |> 
      tidyr::separate(rowname, "_", into = c("assay", "name"), extra = 'merge') |> 
      dplyr::slice(match(tmp$x, name)) %>%
      dplyr::rename(corr_bmi = cor) %>%
      dplyr::rename(p_bmi = p)
    
    tmp <- tmp |> 
      dplyr::left_join(dplyr::select(corr, name, corr_bmi, p_bmi), by = c('x' = 'name'))
    
    hr <- rowAnnotation('Age' = anno_points(tmp$corr_age,
                                            pch = ifelse(tmp$p_age < 0.05, 19, 1),
                                            ylim = c(-1, 1),
                                            gp = gpar(col = ifelse(tmp$corr_age<0, cols[1], cols[4]))),
                        
                        'Bmi' = anno_points(tmp$corr_bmi,
                                            pch = ifelse(tmp$p_bmi < 0.05, 19, 1),
                                            ylim = c(-1, 1),
                                            gp = gpar(col = ifelse(tmp$corr_bmi<0, cols[1], cols[4]))))
  } else {
    hr = NULL
  }
  
  ## Column splits and labels ##---- 
  
  col_split = c(rep("M2", 2), rep("M4", 2), rep("M6", 2))
  
  row_split = ifelse(is.na(tmp$`estimate_visitIdM2:interventionIdK`), 2, 1)
  
  ha = columnAnnotation(estimate = rep(c("time", "modulation MCT"), 3),
                        col = list(estimate = c('time' = cols[3],
                                                'modulation MCT' = cols[8])))
  
  if(!is.null(relabel)){
    
    tmp <- tmp |> 
      # relabel
      dplyr::filter(x %in% relabel$x) |> 
      dplyr::left_join(relabel) |> 
      dplyr::mutate(x = label)
    
  }
  
  ## Draw heatmap ##---- 

  p <- Heatmap(as.matrix(tmp[,c(2, 4, 6, 8, 10, 12)]),
               
               # Row details
               row_labels = tmp$x,
               row_split = row_split,
               #row_title = NULL,
               
               # Row annotation
               right_annotation = hr,
               
               # Column details
               top_annotation = ha,
               column_split = col_split,
               
               # Annotation
               cell_fun = function(j, i, x, y, width, height, fill) {
                 grid.text(as.matrix(tmp[,c(3, 5, 7, 9, 11, 13)])[i, j], x, y, gp = gpar(fontsize = 10))
               },
               row_names_side = 'left',
               column_names_side = 'top',
               # cluster_rows = T,
               row_order = order(tmp$estimate_visitIdM6, decreasing = T),
               cluster_columns = F,
               # column_title = NULL,
               show_row_dend = F, show_column_dend = F,
               show_column_names = F,
               row_title = NULL,
               name = 'estimate\n(scaled)',
               na_col = 'white',
               col = circlize::colorRamp2(breaks = seq(-max(abs(na.omit(tmp[,c(2, 4,6,8,10,12)]))), max(abs(na.omit(tmp[,c(2, 4,6,8,10,12)]))), length.out = 30),
                                          colors = colorRampPalette(c(cols[c(1, 2)], "grey95", cols[c(5, 4)]))(30)
                                          #                            # colors = rev(color("batlow")(5)),
                                          #                            # colors = rev(cols[c(2,5,4, 6, 7)])
               ),
               row_names_gp = grid::gpar(fontsize = 9),
               column_names_gp = grid::gpar(fontsize = 9),
               border_gp = gpar(lwd = 0.5),
               border = T) 
  
  return(p)
  
  
}

