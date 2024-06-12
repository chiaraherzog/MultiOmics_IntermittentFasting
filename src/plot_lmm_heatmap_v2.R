#' @param lmm_data Full model output (lmm helper function)
#' @param cols Colour palette (leave blank for default colours)
#' @param relabel (optional). Data frame supplying labels for the plot. x = variable, label = relabelled variable
#' @param filter_relabel (optional) boolean: if supplying new labels, should only variables provided in that dataframe be plotted? default to TRUE
#' @param colour_assays (optional) if not NULL, list of colours will be used for assays
#' @param methyl_sub default to NULL; filter for methylation data of one tissue only (cervical, buccal, or blood)
#' @param m6_sep Should M6 values be plotted separately to avoid blank spots? default to TRUE
#' @param age_cor should correlation with age be plotted on the side? default to F
#' @return ComplexHeatmap object


plot_lmm_heatmap_v2 <- function(lmm_data,
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
                                relabel_assay = NULL,
                                methyl_sub = NULL,
                                m6_sep = T,
                                age_cor = F,
                                cluster = 'default'){
  
  if(!require("gtools")){
    install.packages("gtools")
  }
  
  tmp <- lmm_data |> 
    dplyr::select(!contains(c("std", "compliancemedium",
                              "bmi_at_consent", 
                              "age_at_consent",
                              "estimate_compliancehigh", "p.value_compliancehigh",
                              "estimate_interventionIdK", "p.value_interventionIdK",
                              "2:interventionIdK", "4:interventionIdK", "6:interventionIdK"))) |> 
    
    # set pval
    dplyr::mutate(across(contains("p.value"), ~ ifelse(is.na(.), "", gtools::stars.pval(.)))) |> 
    dplyr::mutate(across(contains("p.value"), ~ gsub("[.]", "", .)))
  
  # condition: if assay name in x, need to split it for relabelling to work
  if(any(grepl("Composite|Functional|Blood haemogram|Flow cytometry|Body composition", tmp$x))){
    
    tmp <- tmp |> 
      tidyr::separate(x, into = c('assay', "x"),
                      sep = "_", 
                      extra = 'merge')
    
    # methyl subset
    if(!is.null(methyl_sub)){
      tmp <- tmp |> 
        dplyr::filter(grepl(methyl_sub, assay) & grepl("methylation", assay))
    }
    
  }
  
  if(!is.null(relabel)){
    
    if(filter_relabel == T){
      tmp <- tmp |> 
        dplyr::filter(x %in% relabel$x)
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
    
    if(!is.null(relabel_assay)){
      tmp <- tmp |> 
        dplyr::mutate(assay = assay2) |> 
        dplyr::select(-assay2)
    }
    
  }
  
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
  
  row_labs <- ifelse(tmp$`p.value_visitIdM6:compliancehigh` == "*" |tmp$`p.value_visitIdM6:compliancehigh:interventionIdK` == "*" |
                       tmp$`p.value_visitIdM6:compliancehigh` == "**" |tmp$`p.value_visitIdM6:compliancehigh:interventionIdK` == "**" |
                       tmp$`p.value_visitIdM6:compliancehigh` == "***" |tmp$`p.value_visitIdM6:compliancehigh:interventionIdK` == "***" |
                       tmp$`p.value_visitIdM6:compliancehigh` == "." |tmp$`p.value_visitIdM6:compliancehigh:interventionIdK` == ".",
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
                                                                         gp = gpar(col = ifelse(tmp$cor<0, cols[1], cols[4]))))
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
                   grid.text(as.matrix(tmp[,ind_p])[i, j], x, y, gp = gpar(fontsize = 10), just = "centre")
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
  } else {
    p <- Heatmap(as.matrix(tmp[,ind_est]),
                 name = 'estimate\n(scaled)',
                 
                 # Row details
                 row_labels = gt_render(row_labs),
                 row_names_side = 'left',
                 row_split = row_split,
                 cluster_rows = T,
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
                   grid.text(as.matrix(tmp[,ind_p])[i, j], x, y, gp = gpar(fontsize = 10), just = "centre")
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
    
    p <- print(p) + decorate_annotation("Correlation with age\n(baseline)", {
      grid.lines(c(0.5), gp = gpar(col = "grey80",lty = 'dotted'))
    })
    
  }
  
  return(p)
  
}
