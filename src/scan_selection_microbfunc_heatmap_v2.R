#' @param exp which experiment should be plotted
#' @param lmm_data_time Model that provides the time coefficient (no interaction) 
#' @param lmm_data_int Model that provides time*interventionId interaction
#' @param cols Colour palette (leave blank for default colours [recommended], but can be altered if desired)
#' @param selection Data frame supplying annotated Knumbers and pathways selected for plotting.
#' @param pval_threshold Any row with an uncorrected p-value below threshold will be plotted.
#' @return ComplexHeatmap object
#' 
# version made for microbiomes

scan_selection_microbfunc_heatmap_v2 <- function(exp,
                                              lmm_data_time,
                                              lmm_data_int,
                                              cols = c("#1b69a1",
                                                       "#48a0af",
                                                       "#f39668",
                                                       "#ec6669"),
                                              selection,
                                              pval_threshold=1){
  
  ## packages #----
  if(!require("gtools")){
    install.packages("gtools")
  }
  
  require(stringr)
  require(viridisLite)
  
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
  
  ## Step 1: Merge time and MCT estimates, separate assay and feature names and clean up  ##----
  
  tmp1 = lmm_data_time %>%
    filter(str_detect(x, exp)) %>%
    dplyr::select(!contains(c("std.", "age_at_consent", "bmi_at_consent","interventionIdK")))
  
  tmp3 = lmm_data_int %>%
    filter(str_detect(x, exp)) %>%
    dplyr::select(!contains(c("std.", "age_at_consent", "bmi_at_consent",
                              "estimate_interventionIdK",
                              "p.value_interventionIdK"))) %>% 
    dplyr::select(!ends_with(c("visitIdM2", "visitIdM4", "visitIdM6")))
  
  tmp = tmp1 %>% 
    dplyr::full_join(tmp3, by = "x")  %>%
    dplyr::mutate(x = gsub("_clr_","_",x)) %>%
    tidyr::separate(x, into = c('assay', "x"), sep = "_", extra = 'merge')
  
  #tmp$x = gsub(paste0(exp,"_"),"",tmp$x)
  
  ## Step 2: Apply FDR correction, filter out only significant p-values based on p-value threshold and replace p-values by stars and bold significant p-values ##----
  
  tmpfdr = FDRcorr(tmp, append = FALSE)
  
  keep = tmp %>%
    dplyr::filter(rowSums(across(starts_with("p.value"), ~ . < pval_threshold), na.rm = TRUE) > 0) %>%
    dplyr::pull(x)
  
  p_value_columns = grep("^p\\.value", names(tmp), value = TRUE)
  
  for (col in p_value_columns) {
    tmp[[col]] = ifelse(is.na(tmpfdr[[col]]), "",
                        ifelse(tmpfdr[[col]] < 0.05,
                               paste0("<b>", custom_stars_pval(tmp[[col]]), "</b>"), 
                               custom_stars_pval(tmp[[col]])))
  }
  
  tmp = tmp %>% dplyr::filter(x %in% keep)
  
  ## Step 3: Subset features, order and annotate ##----
  
  tmp = tmp %>%
    dplyr::filter(x %in% c(selection$x)) %>% 
    dplyr::left_join(selection) %>%
    dplyr::mutate(x=label)
  
  # order
  tmp = tmp %>% arrange(metabolite.class,Synthesis)
  tmp$metabolite.class = factor(tmp$metabolite.class, levels = unique(tmp$metabolite.class))
  tmp$Synthesis = factor(tmp$Synthesis, levels = unique(tmp$Synthesis))
  
  # Annotation
  cols_class = c("#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","#CC79A7","#FB28BE","#CF83E6","#92B4BA","#88DDD5")
  row_df = tmp %>%
    select(label,Synthesis, metabolite.class) %>%
    column_to_rownames(var="label")
  row_ha = rowAnnotation(df = row_df,
                         col = list(Synthesis = 
                                      setNames(inferno(length(unique(tmp$Synthesis))),unique(tmp$Synthesis)),
                                    metabolite.class = 
                                      setNames(cols_class[1:length(unique(tmp$metabolite.class))],unique(tmp$metabolite.class))),
                         show_annotation_name = c(Synthesis = FALSE,metabolite.class = FALSE))
  
  ## Step 4: Column splits and labels ##---- 
  col_split = rep(c("time", "MCT"), each = 3)
  col_split = factor(col_split, levels = c("time", "MCT"))
  col_labs = rep(c("M2", "M4", "M6"), 2)
  
  # Indices
  tmp = tmp %>% column_to_rownames(var = "x")
  ind_est = grep("estimate", colnames(tmp))  
  ind_p = grep("p.value", colnames(tmp))
  
  ## Step 5: Draw heatmap ##---- 
  
  p = Heatmap(as.matrix(tmp[,ind_est]),
               name = 'estimate\n(scaled)',
               
               # Row details
               row_names_side = 'left',
               cluster_rows = F,
               show_row_dend = F, 
               row_title = NULL,
               
               # Row annotation
               right_annotation = row_ha,
               
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
