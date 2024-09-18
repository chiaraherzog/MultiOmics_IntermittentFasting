#' @param exp which experiment should be plotted
#' @param lmm_data_time Minimal model that provides time coefficient (no interaction) (out_lmm$`Minimal model`)
#' @param lmm_data_compliance Basic model that provides time*compliance interaction (out_lmm$`Basic model with BMI`)
#' @param cols Colour palette (leave blank for default colours [recommended], but can be altered if desired)
#' @param relabel (optional). Data frame supplying labels for the plot. x = variable, label = relabelled variable
#' @param filter_relabel Logical. If `TRUE`, only variables in the relabel data frame will be plotted. Default is `TRUE`.
#' @param age_cor should correlation relative abundance of each feature with age be plotted on the side? default to F
#' @param bmi_cor should correlation relative abundance of each feature with with bmi be plotted on the side? default to F
#' @param buccal_ic_cor should correlation relative abundance of each feature with with buccal_ic be plotted on the side? default to F
#' @param cluster Should rows be clustered? If left blank, will set to 'default' and rows will be ordered by M6 time estimate. If 'cluster' (or other), rows will be clustered.
#' @param pval_threshold Any row with an uncorrected p-value below threshold in the lmm_data_time model will be plotted.
#' @param pval_permutation_filter Logical. If `TRUE`, filters p-values using permutation tests. Default is `TRUE`.
#' @return ComplexHeatmap object
#' 
# version made for microbiomes


# main function
plot_lmm_heatmap_frag_v3 <- function(exp,
                                  lmm_data_time,
                                  lmm_data_compliance,
                                  cols = c("#1b69a1",
                                           "#48a0af",
                                           "#f39668",
                                           "#ec6669"),
                                  relabel = NULL,
                                  filter_relabel = TRUE,
                                  age_cor = F,
                                  bmi_cor = F,
                                  buccal_ic_cor = F,
                                  cluster = 'default',
                                  pval_threshold = 1,
                                  pval_permutation_filter = TRUE){
  ## packages #----
  if(!require("gtools")){
    install.packages("gtools")
  }
  
  ## global helper functions #----
  source("src/FDRcorr.R")
  
  # helper for bolding p-values based on FDR corrected pvalues
  custom_stars_pval <- function(pvals) {
    sapply(pvals, function(pval) {
      if (is.na(pval)) {
        return("")  # Return an empty string for NA values
      } else {
        # Convert p-values to significance stars
        result <- cut(pval, breaks = c(0, 0.001, 0.01, 0.05, 1), 
                      labels = c("***", "**", "*", ""), right = FALSE)
        # Ensure the result is returned as a character
        return(as.character(result))
      }
    })
  }
  
  # helper for heatmap vizualization bold fonts
  cell_fun <- function(j, i, x, y, width, height, fill) {
    full_text <- as.matrix(tmp[, ind_p])[i, j]
    if (grepl("<b>", full_text)) {
      clean_text <- gsub("<b>|</b>", "", full_text)
      grid.text(clean_text, x = x, y = y, 
                gp = gpar(fontsize = 7, fontface = "bold", fontfamily = "Arial Black"), just = "centre")
    } else {
      grid.text(full_text, x = x, y = y, 
                gp = gpar(fontsize = 7, fontface = "plain", col = 'grey20'), just = "centre")
    }
  }
  
  ## Step 1: Filter relevant data frames and p-values, and apply permutation filtering if applicable ##----
  
  tmp1 = lmm_data_time %>%
      filter(str_detect(x, exp)) %>%
      dplyr::select(!contains(c("std.", "age_at_consent", "bmi_at_consent")))
  
  if(pval_permutation_filter) {
    source("src/getPermutationComparisonPval.R")
    var_filter <- tmp1$x[grepl(paste0(relabel$x, collapse = '|'), tmp1$x)]
    
    load("out/outPermMinimal.Rdata")
    tmp1 <- getPermutationComparisonPval(observed = tmp1, random = outPermMinimal, variables = var_filter)
  }
    
  tmp2 = lmm_data_compliance %>%
    filter(str_detect(x, exp)) %>%
    dplyr::select(!contains(c("std", "compliancemedium",
                               "bmi_at_consent", 
                              "age_at_consent",
                               "estimate_compliancehigh", "p.value_compliancehigh"))) %>%
    dplyr::select(!ends_with(c("visitIdM2", "visitIdM4", "visitIdM6")))
    
  if(pval_permutation_filter) {
    load("out/outPerm.Rdata")
    tmp2 <- getPermutationComparisonPval(observed = tmp2, random = outPerm, variables = var_filter)
  }
  
  ## Step 2: Merge time and interaction data, separate assay and feature names, cleanup and filter feature names if applicable ##----
  
  tmp = tmp1 %>% 
    dplyr::full_join(tmp2, by = "x") %>%
    dplyr::mutate(x = gsub("_clr_","_",x)) %>%
    tidyr::separate(x, into = c('assay', "x"), sep = "_", extra = 'merge')
  
  if(!is.null(relabel)) {
    if(filter_relabel) {
      intersect = intersect(tmp$x, relabel$x)
      relabel = relabel[relabel$x %in% intersect,]
      tmp = tmp[match(relabel$x, tmp$x),]
    }}
    
  ## Step 3: Apply FDR correction, filter out only significant p-values among time estimates based on p-value threshold and replace p-values by stars and bold significant p-values ##----
  
  tmpfdr = FDRcorr(tmp, append = FALSE)
  
  keep = tmp %>%
    dplyr::filter(rowSums(across(starts_with("p.value") & ends_with(c("visitIdM2", "visitIdM4", "visitIdM6")), ~ . < pval_threshold), na.rm = TRUE) > 0) %>%
    dplyr::pull(x)
  
  p_value_columns = grep("^p\\.value", names(tmp), value = TRUE)
  
  for (col in p_value_columns) {
    tmp[[col]] = ifelse(is.na(tmpfdr[[col]]), "",
                         ifelse(tmpfdr[[col]] < 0.05,
                                paste0("<b>", custom_stars_pval(tmp[[col]]), "</b>"), 
                                custom_stars_pval(tmp[[col]])))
  }
  
  tmp = tmp %>% dplyr::filter(x %in% keep)
    
  # Step 4: set up indices  ##----
  ind_est <- grep("estimate", colnames(tmp))  
  ind_p <- grep("p.value", colnames(tmp))
  
  ## Step 5: correlation with baseline measurements, if applicable ##----
  
  if(buccal_ic_cor == T) {
    
    load("out/corrIcBaseline_saliva.R")
    corr <- corr |> 
      dplyr::rename(name=rowname) |> 
      dplyr::slice(match(tmp$x, name))
    
    tmp <- tmp |> 
      dplyr::left_join(dplyr::select(corr, name, cor, p), by = c('x' = 'name'))
    
    hr <- rowAnnotation('buccal_ic\n(baseline)' = anno_points(tmp$cor,
                                                              pch = ifelse(tmp$p < 0.05, 19, 1),
                                                              ylim = c(-1, 1),
                                                              gp = gpar(col = ifelse(tmp$cor<0, cols[1], cols[4]))))
    
  } else if (buccal_ic_cor == F & age_cor == T & bmi_cor == F) {
    
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
  } else if (buccal_ic_cor == F & age_cor == F & bmi_cor == T) {
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
  } else if (buccal_ic_cor == F & age_cor == T & bmi_cor == T) {
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
  
  ## Step 6: Column splits and labels ##---- 
  col_split <- factor(rep(c("time", "high\ncompliance"), each = 3), levels = c("time", "high\ncompliance"))
  col_labs <- rep(c("M2", "M4", "M6"), 2)
  
  if(!is.null(relabel)){
    
    tmp <- tmp |> 
      # relabel
      dplyr::filter(x %in% relabel$x) |> 
      dplyr::left_join(relabel) |> 
      dplyr::mutate(x = label)
    
  }
  
  row_labs <- ifelse(grepl("<b>", tmp$`p.value_visitIdM6:compliancehigh`, fixed = TRUE) | 
                       grepl("<b>", tmp$p.value_visitIdM6, fixed = TRUE),
                     paste0("**", tmp$x, "**"),
                     tmp$x)
  
  ## Step 7: Draw heatmap ##---- 

   p <- Heatmap(as.matrix(tmp[,ind_est]),
                name = 'estimate\n(scaled)',
                 
                # Row details
                row_labels = gt_render(row_labs),
                row_names_side = 'left',
                row_order = if(cluster == 'default') order(tmp$estimate_visitIdM6, decreasing = TRUE) else NULL,
                cluster_rows = ifelse(cluster == FALSE, FALSE, TRUE),
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
               cell_fun = cell_fun,
                 
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
  
  return(p)
  
}
