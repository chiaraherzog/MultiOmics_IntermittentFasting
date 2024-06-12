
# Heatmap of features associated with bmi

feat_corr_bmi <- function(dat_long, # features "label"
                                    continous_corr_type = "kendall",
                                    pval = 0.01,
                                    title = "",
                                    change = F){
  
  # correlation bmi
  cor = dat_long %>%
    dplyr::group_by(label) %>%
    dplyr::reframe(cor_bmi = cor(value, bmi_at_consent, method = continous_corr_type),
                   p_bmi = cor.test(value, bmi_at_consent, method = continous_corr_type)$p.value) %>%
    dplyr::ungroup() %>%
    dplyr::filter(p_bmi < pval)
  
  mat <- dat_long %>%
    filter(label %in% cor$label) %>%
    pivot_wider(id_cols = label, names_from = subjectId, values_from = value) %>%
    full_join(cor) %>%
    column_to_rownames(var = "label")
  
  ind_main <- grep("K|I", colnames(mat))
  
  # column annotation
  column_df <- dat_long %>%
    select(subjectId,bmi_at_consent) %>%
    distinct() %>%
    column_to_rownames(var= "subjectId") %>%
    as.data.frame()
  #identical(colnames(mat[,ind_main]),rownames(column_df))
  column_ha = HeatmapAnnotation(df = column_df,
                                col = list(bmi_at_consent = circlize::colorRamp2(
                                             c(round(min(column_df$bmi_at_consent)), 
                                               #round(median(column_df$age_at_consent)), 
                                               (round(max(column_df$bmi_at_consent)) + 1)), 
                                             c("#0FC2C0", "#023535"))
                                ),
                                show_annotation_name = c(bmi_at_consent = FALSE))
  
  # row annotation
  row_ha <- rowAnnotation('BMI' = anno_points(mat$cor_bmi,
                                              pch = ifelse(mat$p_bmi < pval, 19, 1),
                                              ylim = c(-1, 1),
                                              gp = gpar(col = ifelse(mat$cor_bmi<0, cols[1], cols[4]))))
  
  cols = c("#1b69a1","#48a0af", "#f39668","#ec6669")
  
  ind_main <- grep("K|I", colnames(mat))
  
  if (change == F) {
    # set min to 0
    p <- Heatmap(as.matrix(mat[,ind_main]),
                 name = title,
                 
                 # Row details
                 row_names_side = 'left',
                 cluster_rows = T,
                 show_row_dend = F, 
                 row_title = NULL,
                 
                 # Column details
                 cluster_columns = T,
                 show_column_dend = T,
                 
                 # # Annotation
                 top_annotation = column_ha,
                 right_annotation = row_ha,
                 
                 # Colours
                 na_col = 'white',
                 col = circlize::colorRamp2(breaks = seq(0,
                                                         (max(mat[,ind_main], na.rm = TRUE)),
                                                         length.out = 30),
                                            colors = colorRampPalette(c("grey95", cols[c(3, 4)]))(30)
                 ),
                 
                 # Titles
                 row_names_gp = grid::gpar(fontsize = 9),
                 show_column_names = FALSE,
                 
                 # Borders
                 border_gp = gpar(lwd = 0.5),
                 border = T) 
    
    return(p)
  }
}