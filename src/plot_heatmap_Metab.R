#' @description
#' A function to plot (t test) metabolome results plus assay annotation
#' 
#' @param time t test results from high compliance models
#' @param cols Colour palette (leave blank for default colours [recommended], but can be altered if desired)
#' @param relabel (optional). Data frame supplying labels for the plot. x = rowname, label = relabelled rowname
#' @param filter_relabel (optional) boolean: if supplying new labels, should only variables provided in that dataframe be plotted? default to TRUE
#' @param relabel_assay (optional). Set to TRUE if you would like to rename the assay. Assay names should be provided in the relabel parameter, under column name 'assay2'
#' @param colour_assays (optional) if not NULL, list of default colours will be used for assays. This only works for clinical variables so far.
#' @param cluster Should rows be clustered? If left blank, will set to 'default' and rows will be ordered by M6 time estimate. If 'cluster' (or other), rows will be clustered.
#' @param change_baseline Should the change be converted in % from baseline mean?
#' @return ComplexHeatmap object

plot_heatmap_Metab <- function(time,
                                cols = c("#1b69a1",
                                         "#48a0af",
                                         "#f39668",
                                         "#ec6669"),
                                relabel = NULL,
                                filter_relabel = T,
                                colour_assays = list(class = c("Alkaloids" = "#1B69A1",
                                                               "Benzenoids" = "#3D93AB",
                                                               "Carbohydrates" = "#5EABAB",
                                                               "Fatty Acyls" = "#9CAA93",
                                                               "Nucleic acids" = "#F29068",
                                                               "Organic acids" = "#EC6B68",
                                                               "Organic nitrogen compounds" = "#CC6476",
                                                               "Organic oxygen compounds" = "#A34B8A",
                                                               "Organoheterocyclic compounds" = "#7A3B9D",
                                                               "Polyketides" = "#5F70A8")),
                                relabel_assay = T,
                                cluster = 'default',
                               change_baseline = F,
                               padj = F){
  
  if(!require("gtools")){
    install.packages("gtools")
  }
  
  if(change_baseline == T){
    time <- time |> 
      dplyr::mutate(across(c(M2, M4, M6), ~ . / mean_baseline*100))
  }
  
  # Filter metabolites
  tmp1 <- time |> 
    dplyr::select(assay, rowname, M2, p_M2, M4, p_M4, M6, p_M6)
  
  tmp <- tmp1 |> 
    dplyr::filter(rowname %in% relabel$x) |> 
    dplyr::left_join(relabel, by = c("rowname" = "x",
                                     'assay' = 'assay')) |> 
    dplyr::filter(! p_M2 %in% c("", " ") |! p_M4 %in% c("", " ") | ! p_M6 %in% c("", " ")) |> 
    dplyr::arrange(assay2)
  
  
  if(padj == T){
    tmp <- tmp |> dplyr::mutate(across(starts_with("p_") | starts_with("diff_p"), ~ p.adjust(., method = 'fdr')))
  }
  
  tmp <- tmp |>  
    dplyr::mutate(across(starts_with("p_") | starts_with("diff_p"), ~ ifelse(is.na(.), "", 
                                                                             ifelse(. < 0.05,
                                                                                    paste0("<b>",gtools::stars.pval(.),"</b>"), 
                                                                                    "")))) |>
    dplyr::mutate(across(starts_with("p_") | starts_with("diff_p"), ~ gsub("<b>[.]</b>", "", .)))
      
  
  # Indices
  ind_est <- grep("^M2|^M4|^M6|^diff_M2|^diff_M4|^diff_M6", colnames(tmp))  
  ind_p <- grep("^p", colnames(tmp))
  row_labs <- ifelse(grepl("<b>", tmp$p_M6, fixed = TRUE),
                     paste0("**", tmp$label, "**"),
                     tmp$label)
  
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
  
  ha <- rowAnnotation(class = tmp$assay2,
                      col = colour_assays)
  name = ifelse(change_baseline == T, 'change\n(% from\nbaseline)', 'change\n(absolute)')
  
  if(cluster == F){
    
    p <- Heatmap(as.matrix(tmp[,ind_est]),
                 name = name,
                 
                 # Row details
                 row_labels = gt_render(row_labs),
                 row_names_side = 'left',
                 cluster_rows = F,
                 # row_split = row_split,
                 # row_order = order(tmp$estimate_visitIdM6, decreasing = T),
                 show_row_dend = F, 
                 row_title = NULL,
                 
                 # Row annotation
                 left_annotation = ha,
                 
                 # Column details
                 cluster_columns = F,
                 show_column_dend = F,
                 
                 # Row side colours
                 
                 
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
    
  } else if (cluster == T){
    p <- Heatmap(as.matrix(tmp[,ind_est]),
                 name = name,
                 
                 # Row details
                 row_labels = gt_render(row_labs),
                 row_names_side = 'left',
                 cluster_rows = T,
                 # row_split = row_split,
                 # row_order = order(tmp$estimate_visitIdM6, decreasing = T),
                 show_row_dend = F, 
                 row_title = NULL,
                 
                 # Row annotation
                 left_annotation = ha,
                 
                 # Column details
                 cluster_columns = F,
                 show_column_dend = F,
                 
                 # Row side colours
                 
                 
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
  } else {
    p <- Heatmap(as.matrix(tmp[,ind_est]),
                 name = name,
                 
                 # Row details
                 row_labels = gt_render(row_labs),
                 row_names_side = 'left',
                 cluster_rows = F,
                 # row_split = row_split,
                 row_order = order(tmp$M6, decreasing = T),
                 show_row_dend = F, 
                 row_title = NULL,
                 
                 # Row annotation
                 left_annotation = ha,
                 
                 # Column details
                 cluster_columns = F,
                 show_column_dend = F,
                 
                 # Row side colours
                 
                 
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
  }
    

  return(p)
  
}
