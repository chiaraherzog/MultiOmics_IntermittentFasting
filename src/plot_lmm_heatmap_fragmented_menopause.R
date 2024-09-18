#' @param exp which experiment should be plotted
#' @param lmm_data_overall General model with menopauze term (no interaction) (out_lmm$`Menopause general`)
#' @param lmm_data_int Data from model with additional interaction term. (out_lmm$`Menopause interaction`)
#' @param cols Colour palette (leave blank for default colours [recommended], but can be altered if desired)
#' @param FDR should p-values be corrected? default to F
#' @return ComplexHeatmap object
#' 
# version made for additional checks on menopauzal status
# models run for high compliance individuals only
# option to FDR correct p-values, for each term (and category) separately 
# all terms with p.corrected < 0.05 are plotted.
# for uncorrected p-values, all rows with p < 0.05 are plotted

# main function
plot_lmm_heatmap_frag <- function(exp,
                                  lmm_data_overall,
                                  lmm_data_int,
                                  cols = c("#1b69a1",
                                           "#48a0af",
                                           "#f39668",
                                           "#ec6669"),
                                  FDR = T){
  
  if(!require("gtools")){
    install.packages("gtools")
  }
  
  library(stringr)
  library(purrr)
  library(openxlsx)
  library(rcartocolor)
    
  tmp1 <- lmm_data_overall |> 
    dplyr::filter(purrr::map_lgl(x, ~ any(stringr::str_detect(.x, exp)))) |> 
    dplyr::select(x|contains(c("mpstatrsyes"))) |> 
    dplyr::select(!contains(c("std.")))
    
  tmp3 <- lmm_data_int |> 
    dplyr::filter(purrr::map_lgl(x, ~ any(stringr::str_detect(.x, exp)))) %>%
    dplyr::select(!contains(c("std.", "age_at_consent", "bmi_at_consent","interventionIdK",
                              "estimate_mpstatrsyes",
                              "p.value_mpstatrsyes"))) %>% 
    dplyr::select(!ends_with(c("visitIdM2", "visitIdM4", "visitIdM6", "visitIdM2_adj",
                               "visitIdM4_adj", "visitIdM6_adj",
                               "group")))
    
  tmp <- tmp1 |> 
    dplyr::full_join(tmp3, by = "x") |> 
    dplyr::filter(!grepl("families_clr|ASVs_ASV|Saliva nuclear", x))
    
  
  p_value_columns <- grep("^p\\.value", names(tmp), value = TRUE)
  
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
  
  if(fdr){
    tmpfdr <- FDRcorr(tmp, append = F)
    
    for (col in p_value_columns) {
      tmp[[col]] <- ifelse(is.na(tmpfdr[[col]]), "",
                           ifelse(tmpfdr[[col]] < 0.05,
                                  paste0("<b>", custom_stars_pval(tmp[[col]]), "</b>"), 
                                  custom_stars_pval(tmp[[col]])))
    }
    
  } else {
    
    for (col in p_value_columns) {
      tmp[[col]] <- ifelse(is.na(tmp[[col]]), "",
                           ifelse(tmp[[col]] < 0.05,
                                  paste0("<b>", custom_stars_pval(tmp[[col]]), "</b>"), 
                                  custom_stars_pval(tmp[[col]])))
    }
    
  }
  

  # Filter only those where at least one is significant
  tmp <- tmp |> 
    filter_at(vars(starts_with("p.value")), any_vars(grepl("**", ., fixed = TRUE)))
    
  ## Split feature names into label/ome, fix names, order by ome ##----
  tmp <- tmp |> 
    dplyr::mutate(x = gsub("_clr", "", x)) |>
    tidyr::separate(x, "_", into = c("assay", "x"), extra = 'merge')
  
  load("src/vars.Rdata")
  
  tmp <- tmp |> dplyr::left_join(dplyr::select(vars, x, label, assay, assay2)) |> 
    dplyr::mutate(x = label,
                  
                  assay2 = case_when(assay2 %in% c('Blood test') | assay == 'Blood haemogram' ~ "Routine bloods",
                                                   assay2 %in% c("Functional exercise capacity",
                                                                 "Body weight and composition") | assay %in% c("Body composition", "Functional sports exam") ~ "Functional clinical features",
                                                   grepl("Flow cytometry", assay) ~ "Immune",
                                                   grepl("cervical", assay) ~ "Cervical methylation",
                                                   grepl("buccal", assay) ~ "Buccal methylation",
                                                   grepl("blood", assay) & grepl("methylation", assay) ~ "Blood methylation",
                                                   grepl("ASV", assay) & grepl("Saliva", assay) ~ "Saliva microbiome (ASV)",
                                                   grepl("families", assay) & grepl("Saliva", assay) ~ "Saliva microbiome (family)",
                                                   grepl("ASV", assay) & grepl("Stool", assay) ~ "Stool microbiome (ASV)",
                                                   grepl("families", assay) & grepl("Stool", assay) ~ "Stool microbiome (family)",
                                                   grepl("magnetic", assay) & grepl("Urine", assay) ~ "Urine metabolome",
                                     grepl("magnetic", assay) & grepl("Saliva", assay) ~ "Saliva metabolome",
                                                   TRUE ~ NA),
  
                  assay = assay2) |> 
    dplyr::select(-c(assay2, label)) |> 
    dplyr::relocate(assay, x)
  
  # order by ome
  tmp <- tmp |> dplyr::mutate(assay = factor(assay, levels = c("Functional clinical features",
                                                               "Blood methylation",
                                                               "Buccal methylation",
                                                               "Cervical methylation",
                                                               "Immune",
                                                               "Saliva microbiome (family)",
                                                               "Stool microbiome (family)",
                                                               "Urine metabolome"))) |> 
    droplevels() |> 
    arrange(assay)
  
  xcol = list(assay = c("Functional clinical features" = "#1b69a1",
                        "Blood methylation" = "#48a0af",
                        "Buccal methylation" = "#71b5a9",
                        "Cervical methylation" = "#f39668",
                        "Immune" = "#ec6669",
                        "Saliva microbiome (family)" = "#bd647d",
                        "Stool microbiome (family)" = "#832c9b",
                        "Urine metabolome" =  "#5f70a8"))
  
  ## Column splits and labels ##---- 
  col_split <- c("menopause", "menopause:time","menopause:time","menopause:time")
  col_split <- factor(col_split, levels = c("menopause", "menopause:time"))
  col_labs <- c("","M2", "M4", "M6")
  
  row_labs <- ifelse(tmp$`p.value_visitIdM6:mpstatrsyes` == "*" |
                     tmp$`p.value_visitIdM6:mpstatrsyes` == "**" |
                     tmp$`p.value_visitIdM6:mpstatrsyes` == "***" |
                     tmp$`p.value_visitIdM6:mpstatrsyes` == ".",
                     paste0("**", tmp$x, "**"),
                     tmp$x)
  
  row_labs <- ifelse(grepl("<b>", tmp$`p.value_visitIdM6:mpstatrsyes`, fixed = TRUE) | 
                       grepl("<b>", tmp$`p.value_visitIdM4:mpstatrsyes`, fixed = TRUE) |
                       grepl("<b>", tmp$`p.value_visitIdM2:mpstatrsyes`, fixed = TRUE) |
                       grepl("<b>", tmp$p.value_mpstatrsyes, fixed = TRUE),
                     paste0("**", tmp$x, "**"),  # FDR significant (bold)
                            tmp$x)
  
  # ome annotation
  ha <- rowAnnotation(assay = tmp$assay,
                      col = xcol)

  ## Draw heatmap ##---- 
  
  # Indices
  ind_est <- grep("estimate", colnames(tmp))  
  ind_p <- grep("p.value", colnames(tmp))
  
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

  p <- Heatmap(as.matrix(tmp[,ind_est]),
               name = 'estimate\n(scaled)',
                
               # Row details
               row_labels = gt_render(row_labs),
               row_names_side = 'left',
                cluster_rows = F,
                cluster_row_slices = F,
               show_row_dend = F, 
               row_title = NULL,
               
                # Row annotation
               left_annotation = ha,
               
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
  
  
  return(p)
  
}
