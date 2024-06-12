#' @param lmm_data_time Minimal model that provides time coefficient (no interaction) (out_lmm$`Minimal model`)
#' @param lmm_data_compliance Basic model that provides time*compliance interaction (out_lmm$`Basic model with BMI`)
#' @param lmm_data_int High compliance intervention output. (out_lmm$`Intervention (higher compliance)`)
#' @param cols Colour palette (leave blank for default colours [recommended], but can be altered if desired)
#' @param relabel (optional). Data frame supplying labels for the plot. x = variable, label = relabelled variable
#' @param filter_relabel (optional) boolean: if supplying new labels, should only variables provided in that dataframe be plotted? default to TRUE
#' @param relabel_assay (optional). Set to TRUE if you would like to rename the assay. Assay names should be provided in the relabel parameter, under column name 'assay2'
#' @param colour_assays (optional) if not NULL, list of default colours will be used for assays. This only works for clinical variables so far.
#' @param methyl_sub default to NULL; filter for methylation data of one tissue only (cervical, buccal, or blood)
#' @param m6_sep Should M6 values be plotted separately to avoid blank spots? default to TRUE
#' @param age_cor should correlation with age be plotted on the side? default to F
#' @param mark_age_cor description
#' @param cluster Should rows be clustered? If left blank, will set to 'default' and rows will be ordered by M6 time estimate. If 'cluster' (or other), rows will be clustered.
#' @return ComplexHeatmap object


plot_lmm_heatmap_v3 <- function(lmm_data_time,
                                lmm_data_compliance,
                                lmm_data_int,
                                cols = c("#1b69a1",
                                         "#48a0af",
                                         "#f39668",
                                         "#ec6669"),
                                relabel = NULL,
                                filter_relabel = T,
                                colour_assays = list(assay = c("Blood haemogram" = "#ec6669",
                                                               "Body weight and composition" = "#832c9b",
                                                               "Spirometry" = "#bd647d",
                                                               "Functional exercise capacity" = "#f39668",
                                                               "Blood test" = '#71b5a9',
                                                               "Skin histology and TEWL" = "#5f70a8",
                                                               "Vascular features" = '#0d49a1')),
                                relabel_assay = F,
                                methyl_sub = NULL,
                                m6_sep = T,
                                age_cor = F,
                                mark_age_cor = F,
                                cluster = 'default',
                                padj = F){
  
  if(!require("gtools")){
    install.packages("gtools")
  }
  
  tmp1 <- lmm_data_time |> 
    dplyr::select(!contains(c("std.", "age_at_consent", "bmi_at_consent")))
  
  tmp2 <- lmm_data_compliance |> 
    dplyr::select(!contains(c("std", "compliancemedium",
                              "bmi_at_consent", 
                              "age_at_consent",
                              "estimate_compliancehigh", "p.value_compliancehigh"))) |> 
    dplyr::select(!ends_with(c("visitIdM2", "visitIdM4", "visitIdM6")))
  
  tmp3 <- lmm_data_int |> 
    dplyr::select(!contains(c("std.", "age_at_consent", "bmi_at_consent", "compliancemedium",
                              "estimate_interventionIdK",
                              "p.value_interventionIdK"))) |> 
    dplyr::select(!ends_with(c("visitIdM2", "visitIdM4", "visitIdM6")))
  
  tmp <- tmp1 |> 
    dplyr::left_join(tmp2) |> 
    dplyr::left_join(tmp3) |> 
      tidyr::separate(x, into = c('assay', "x"),
                      sep = "_", 
                      extra = 'merge')
  

  # methyl subset
  if(!is.null(methyl_sub)){
      tmp <- tmp |> 
        dplyr::filter(grepl(methyl_sub, assay) & grepl("methylation", assay))
    }
  
  if(!is.null(relabel)){
    
    if(filter_relabel == T){
      intersect <- intersect(tmp$x, relabel$x)
      relabel <- relabel[relabel$x %in% intersect,]
      tmp <- tmp[match(relabel$x, tmp$x),]
    }
    
    if(age_cor == T){
      load("out/corrAgeBaseline.R")
      corr <- corr |> 
        tidyr::separate(rowname, "_", into = c("assay", "name"), extra = 'merge') |> 
        dplyr::slice(match(tmp$x, name))
      
      tmp <- tmp |> 
        dplyr::left_join(dplyr::select(corr, name, cor, p), by = c('x' = 'name'))
    }
    
    tmp <- tmp |> 
      dplyr::left_join(relabel) |> 
      dplyr::mutate(x = label) |> 
      dplyr::select(-label)
    
    if(relabel_assay == T){
      tmp <- tmp |> 
        dplyr::mutate(assay = assay2) |> 
        dplyr::select(-assay2)
    }
    
  }
  
  if(padj == T){
    tmp <- tmp |> dplyr::mutate(across(contains("p.value"), ~ p.adjust(., method = 'fdr')))
  }
  
  tmp <- tmp |> 
    # set pval
    dplyr::mutate(across(contains("p.value"), ~ ifelse(is.na(.), "", gtools::stars.pval(.)))) |> 
    dplyr::mutate(across(contains("p.value"), ~ gsub("[.]", "", .)))
  
  
  # Indices
  ind_est <- grep("estimate", colnames(tmp))  
  ind_p <- grep("p.value", colnames(tmp))
  
  # Column splits and labels
  col_split <- rep(c("time", "high\ncompliance", "MCT"), each = 3)
  col_split <- factor(col_split, levels = c("time", "high\ncompliance", "MCT"))
  col_labs <- rep(c("M2", "M4", "M6"), 3)
  
  
  # Row splits and labels
  if(m6_sep == T){
    row_split = ifelse(is.na(tmp$`estimate_visitIdM2`), 2, 1)
  } else {
    row_split = rep(1, nrow(tmp))
  }
  
  row_labs <- ifelse(tmp$`p.value_visitIdM6:compliancehigh` == "*" |tmp$`p.value_visitIdM6:interventionIdK` == "*" |
                       tmp$`p.value_visitIdM6:compliancehigh` == "**" |tmp$`p.value_visitIdM6:interventionIdK` == "**" |
                       tmp$`p.value_visitIdM6:compliancehigh` == "***" |tmp$`p.value_visitIdM6:interventionIdK` == "***" |
                       tmp$`p.value_visitIdM6:compliancehigh` == "." |tmp$`p.value_visitIdM6:interventionIdK` == ".",
                     paste0("**", tmp$x, "**"),
                     tmp$x)
  
  if(is.list(colour_assays)){
    ha <- rowAnnotation(assay = tmp$assay,
                        col = colour_assays)
  } else {
    ha <- NULL
  }
  
  if(age_cor == T){
    hr <- rowAnnotation('Correlation with age\n(baseline)' = anno_points(tmp$cor,
                                                                         pch = ifelse(tmp$p < 0.05, 19, 1),
                                                                         ylim = c(-1, 1),
                                                                         gp = gpar(col = ifelse(tmp$cor<0, cols[1], cols[4]))),
                        
                        'Opposite change\nwith\nintervention' = anno_text(ifelse(sign(tmp$cor) != sign(tmp$estimate_visitIdM6) & tmp$p.value_visitIdM6!=" ",
                                                                                 "←", "")))
  } else {
    hr = NULL
  }
  
  if(cluster == 'default'){
    p <- Heatmap(as.matrix(tmp[,ind_est]),
                 name = 'estimate\n(scaled)',
                 
                 # Row details
                 row_labels = gt_render(row_labs),
                 row_names_side = 'left',
                 row_split = row_split,
                 row_order = order(tmp$estimate_visitIdM6, decreasing = T),
                 show_row_dend = F, 
                 row_title = NULL,
                 
                 # Row annotation
                 left_annotation = ha,
                 right_annotation = hr,
                 
                 # Column details
                 column_labels = col_labs,
                 column_title_gp = grid::gpar(fontsize = 10),
                 column_split = col_split,
                 cluster_columns = F,
                 show_column_dend = F,
                 
                 # Row side colours
                 
                 
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
  } else if(cluster==F){
    
    p <- Heatmap(as.matrix(tmp[,ind_est]),
                 name = 'estimate\n(scaled)',
                 
                 # Row details
                 row_labels = gt_render(row_labs),
                 row_names_side = 'left',
                 row_split = row_split,
                 cluster_rows = F,
                 cluster_row_slices = F,
                 show_row_dend = F, 
                 row_title = NULL,
                 
                 # Row annotation
                 left_annotation = ha,
                 right_annotation = hr,
                 
                 # Column details
                 column_labels = col_labs,
                 column_title_gp = grid::gpar(fontsize = 10),
                 column_split = col_split,
                 cluster_columns = F,
                 show_column_dend = F,
                 
                 # Row side colours
                 
                 
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
    if(age_cor == T){
      p <- print(p) + decorate_annotation("Correlation with age\n(baseline)", {
        grid.lines(c(0.5), gp = gpar(col = "grey80",lty = 'dotted'))
      })
    }
    
  } else {
    p <- Heatmap(as.matrix(tmp[,ind_est]),
                 name = 'estimate\n(scaled)',
                 
                 # Row details
                 row_labels = gt_render(row_labs),
                 row_names_side = 'left',
                 row_split = row_split,
                 cluster_rows = T,
                 cluster_row_slices = F,
                 show_row_dend = F, 
                 row_title = NULL,
                 
                 # Row annotation
                 left_annotation = ha,
                 right_annotation = hr,
                 
                 # Column details
                 column_labels = col_labs,
                 column_title_gp = grid::gpar(fontsize = 10),
                 column_split = col_split,
                 cluster_columns = F,
                 show_column_dend = F,
                 
                 # Row side colours
                 
                 
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
    if(age_cor == T){
    p <- print(p) + decorate_annotation("Correlation with age\n(baseline)", {
      grid.lines(c(0.5), gp = gpar(col = "grey80",lty = 'dotted'))
    })
    }
    
  }
  
  return(p)
  
}
