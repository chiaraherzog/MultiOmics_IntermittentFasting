#' @param exp which experiment should be plotted
#' @param lmm_data_time Minimal model that provides time coefficient (no interaction) (out_lmm$`Minimal model`)
#' @param lmm_data_compliance Basic model that provides time*compliance interaction (out_lmm$`Basic model with BMI`)
#' @param lmm_data_int High compliance intervention output. (out_lmm$`Intervention (higher compliance)`)
#' @param cols Colour palette (leave blank for default colours [recommended], but can be altered if desired)
#' @param relabel (optional). Data frame supplying labels for the plot. x = variable, label = relabelled variable
#' @param age_cor should correlation with age be plotted on the side? default to F
#' @param bmi_cor should correlation with bmi be plotted on the side? default to F
#' @param FDR should p-values be corrected? default to F
#' @param cluster Should rows be clustered? If left blank, will set to 'default' and rows will be ordered by M6 time estimate. If 'cluster' (or other), rows will be clustered.
#' @return ComplexHeatmap object
#' 
# version made for microbiomes
# option to plot age and/or bmi correlations at baseline
# option to FDR correct p-values, for each term (and category) separately 
# all terms with p.corrected < 0.05 are plotted.
# for uncorrected p-values, all rows with p < 0.01 are plotted

# main function
plot_lmm_heatmap_frag <- function(exp,
                                  lmm_data_time,
                                  lmm_data_compliance,
                                  lmm_data_int,
                                  cols = c("#1b69a1",
                                           "#48a0af",
                                           "#f39668",
                                           "#ec6669"),
                                  relabel = NULL,
                                  age_cor = F,
                                  bmi_cor = F,
                                  FDR = F,
                                  cluster = 'default'){
  
  if(!require("gtools")){
    install.packages("gtools")
  }
  
  source("src/FDR_correct.R")
  
  ## load output lmm models and correct p-values if requested, keep only significant results ##----
  
  if (FDR == F) {
    
    tmp1 <- lmm_data_time %>%
      filter(str_detect(x, exp)) %>%
      dplyr::select(!contains(c("std.", "age_at_consent", "bmi_at_consent")))
    
    tmp2 <- lmm_data_compliance %>%
      filter(str_detect(x, exp)) %>%
      dplyr::select(!contains(c("std", "compliancemedium",
                                "bmi_at_consent", 
                                "age_at_consent",
                                "estimate_compliancehigh", "p.value_compliancehigh"))) %>%
      dplyr::select(!ends_with(c("visitIdM2", "visitIdM4", "visitIdM6")))
    
    tmp3 <- lmm_data_int %>%
      filter(str_detect(x, exp)) %>%
      dplyr::select(!contains(c("std.", "age_at_consent", "bmi_at_consent", "compliancemedium",
                                "estimate_interventionIdK",
                                "p.value_interventionIdK"))) %>% 
      dplyr::select(!ends_with(c("visitIdM2", "visitIdM4", "visitIdM6")))
    
    tmp <- tmp1 %>% 
      dplyr::full_join(tmp2, by = "x") %>% 
      dplyr::full_join(tmp3, by = "x")
    
    tmp$x = gsub(paste0(exp,"_"),"",tmp$x)
    
    # filter out only significant rows and replace pval with * symbols
    tmp <- tmp %>%
      filter(rowSums(across(starts_with("p.value"), ~ . < 0.01), na.rm = TRUE) > 0) %>% 
      dplyr::mutate(across(contains("p.value"), ~ ifelse(is.na(.), "", gtools::stars.pval(.)))) |> 
      dplyr::mutate(across(contains("p.value"), ~ gsub("[.]", "", .)))

  } else if (FDR == T){
    
    # simple model FDR correct
    tmp1 <- lmm_data_time %>%
      filter(str_detect(x, exp))
    
    tmp1_p <- tmp1 %>%
      dplyr::select(x,contains("p.value")) %>%
      pivot_longer(cols = contains("p.value"),names_to = "estimate", values_to = "p_value") %>%
      FDR_correct() %>%
      pivot_wider(id_cols = x, names_from = estimate, values_from = p_value)
    
    tmp1 <- tmp1 %>%
      select(!contains("p.value")) %>%
      left_join(tmp1_p, by = "x") %>%
      dplyr::select(!contains(c("std.", "age_at_consent", "bmi_at_consent")))
    
    # compliance model FDR correct
    tmp2 <- lmm_data_compliance %>%
      filter(str_detect(x, exp))
    
    tmp2_p <- tmp2 %>%
      dplyr::select(x,contains("p.value")) %>%
      pivot_longer(cols = contains("p.value"),names_to = "estimate", values_to = "p_value") %>%
      FDR_correct() %>%
      pivot_wider(id_cols = x, names_from = estimate, values_from = p_value)
    
    tmp2 <- tmp2 %>%
      select(!contains("p.value")) %>%
      left_join(tmp2_p, by = "x") %>%
      dplyr::select(!contains(c("std", "compliancemedium",
                                "bmi_at_consent", 
                                "age_at_consent",
                                "estimate_compliancehigh", "p.value_compliancehigh"))) %>%
      dplyr::select(!ends_with(c("visitIdM2", "visitIdM4", "visitIdM6")))
    
    # MCT high compliance part FDR correct
    tmp3 <- lmm_data_int %>%
      filter(str_detect(x, exp))
    
    tmp3_p <- tmp3 %>%
      dplyr::select(x,contains("p.value")) %>%
      pivot_longer(cols = contains("p.value"),names_to = "estimate", values_to = "p_value") %>%
      FDR_correct() %>%
      pivot_wider(id_cols = x, names_from = estimate, values_from = p_value)
    
    tmp3 <- tmp3 %>%
      select(!contains("p.value")) %>%
      left_join(tmp3_p, by = "x") %>%
      dplyr::select(!contains(c("std.", "age_at_consent", "bmi_at_consent", "compliancemedium",
                                "estimate_interventionIdK",
                                "p.value_interventionIdK"))) %>% 
      dplyr::select(!ends_with(c("visitIdM2", "visitIdM4", "visitIdM6")))
    
    # combine
    tmp <- tmp1 %>% 
      dplyr::full_join(tmp2, by = "x") %>% 
      dplyr::full_join(tmp3, by = "x")
    
    tmp$x = gsub(paste0(exp,"_"),"",tmp$x)
    
    # filter out only significant rows and replace pval with * symbols
    tmp <- tmp %>%
      filter(rowSums(across(starts_with("p.value"), ~ . < 0.05)) > 0, na.rm = TRUE) %>% 
      dplyr::mutate(across(contains("p.value"), ~ ifelse(is.na(.), "", gtools::stars.pval(.)))) |> 
      dplyr::mutate(across(contains("p.value"), ~ gsub("[.]", "", .)))
  }
  
  # Indices
  ind_est <- grep("estimate", colnames(tmp))  
  ind_p <- grep("p.value", colnames(tmp))
  
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
  col_split <- rep(c("time", "high\ncompliance", "MCT"), each = 3)
  col_split <- factor(col_split, levels = c("time", "high\ncompliance", "MCT"))
  col_labs <- rep(c("M2", "M4", "M6"), 3)
  
  if(!is.null(relabel)){
    
    tmp <- tmp |> 
      # relabel
      dplyr::filter(x %in% relabel$x) |> 
      dplyr::left_join(relabel) |> 
      dplyr::mutate(x = label)
    
  }
  
  row_labs <- ifelse(tmp$`p.value_visitIdM6:compliancehigh` == "*" |tmp$`p.value_visitIdM6:interventionIdK` == "*" |
                       tmp$`p.value_visitIdM6:compliancehigh` == "**" |tmp$`p.value_visitIdM6:interventionIdK` == "**" |
                       tmp$`p.value_visitIdM6:compliancehigh` == "***" |tmp$`p.value_visitIdM6:interventionIdK` == "***" |
                       tmp$`p.value_visitIdM6:compliancehigh` == "." |tmp$`p.value_visitIdM6:interventionIdK` == ".",
                     paste0("**", tmp$x, "**"),
                     tmp$x)
  
  ## Draw heatmap ##---- 
  
  if(cluster == 'default'){
    p <- Heatmap(as.matrix(tmp[,ind_est]),
                 name = 'estimate\n(scaled)',
                 
                 # Row details
                 row_labels = gt_render(row_labs),
                 row_names_side = 'left',
                 row_order = order(tmp$estimate_visitIdM6, decreasing = T),
                 show_row_dend = F, 
                 row_title = NULL,
                 
                 # Row annotation
                 right_annotation = hr,
                 
                 # Column details
                 column_labels = col_labs,
                 column_title_gp = grid::gpar(fontsize = 10),
                 column_split = col_split,
                 cluster_columns = F,
                 show_column_dend = F,
                 
                 # Annotation
                 cell_fun = function(j, i, x, y, width, height, fill) {
                   grid.text(as.matrix(tmp[,ind_p])[i, j], x, y, gp = gpar(fontsize = 7), just = "centre")
                 },
                 
                 # Colours
                 na_col = 'white',
                 col = circlize::colorRamp2(breaks = seq(-max(abs(na.omit(tmp[,ind_est]))),
                                                         max(abs(na.omit(tmp[,ind_est]))),
                                                         length.out = 30),
                                            colors = colorRampPalette(c(cols[c(1, 2)],
                                                                        "grey95", cols[c(3, 4)]))(30)),
                 
                 # Titles
                 row_names_gp = grid::gpar(fontsize = 9),
                 column_names_gp = grid::gpar(fontsize = 9),
                 
                 # Borders
                 border_gp = gpar(lwd = 0.5),
                 border = T)  
  } else {
    p <- Heatmap(as.matrix(tmp[,ind_est]),
                 name = 'estimate\n(scaled)',
                 
                 # Row details
                 row_labels = gt_render(row_labs),
                 row_names_side = 'left',
                 cluster_rows = T,
                 cluster_row_slices = F,
                 show_row_dend = F, 
                 row_title = NULL,
                 
                 # Row annotation
                 right_annotation = hr,
                 
                 # Column details
                 column_labels = col_labs,
                 column_title_gp = grid::gpar(fontsize = 10),
                 column_split = col_split,
                 cluster_columns = F,
                 show_column_dend = F,
                 
                 # Annotation
                 cell_fun = function(j, i, x, y, width, height, fill) {
                   grid.text(as.matrix(tmp[,ind_p])[i, j], x, y, gp = gpar(fontsize = 7), just = "centre")
                 },
                 
                 # Colours
                 na_col = 'white',
                 col = circlize::colorRamp2(breaks = seq(-max(abs(na.omit(tmp[,ind_est]))),
                                                         max(abs(na.omit(tmp[,ind_est]))),
                                                         length.out = 30),
                                            colors = colorRampPalette(c(cols[c(1, 2)],
                                                                        "grey95", cols[c(3, 4)]))(30)
                                            #                            # colors = rev(color("batlow")(5)),
                                            #                            # colors = rev(cols[c(2,5,4, 6, 7)])
                 ),
                 
                 # Titles
                 row_names_gp = grid::gpar(fontsize = 9),
                 column_names_gp = grid::gpar(fontsize = 9),
                 
                 # Borders
                 border_gp = gpar(lwd = 0.5),
                 border = T)  
  }
  
  return(p)
  
}
