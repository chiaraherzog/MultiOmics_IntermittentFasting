plot_lmm_heatmap_baseline_immune <- function(lmm_data,
                                      cols = cols,
                                      type = 'compliance',
                                      populations = NULL,
                                      corr_baseline = NULL){
  
  # Relabel
  lmm_data <- lmm_data |> 
    dplyr::mutate(name = x) |> 
    dplyr::left_join(dplyr::select(populations, name, pop_name, type)) |> 
    dplyr::left_join(dplyr::select(corr_baseline, name, cor, p.adj))
    
  
  if(type == 'compliance'){
    
    tmp <- lmm_data |> 
      dplyr::select(name, 
                    `estimate_visitIdM2:compliancehigh`, `p.value_visitIdM2:compliancehigh`,
                    `estimate_visitIdM4:compliancehigh`, `p.value_visitIdM4:compliancehigh`,
                    `estimate_visitIdM6:compliancehigh`, `p.value_visitIdM6:compliancehigh`,
                    cor, p.adj,
                    pop_name, type) |> 
      dplyr::mutate(across(contains("p.value"), ~ ifelse(.<0.05, "*", ""))) |> 
      dplyr::rename_at(vars(contains('estimate_')), ~ gsub("estimate_visitId|:compliancehigh", "", .)) |> 
      dplyr::filter(!is.na(pop_name))
    
    if(!is.null(corr_baseline)){
      # annotation 
      left_anno <- rowAnnotation('Correlation with age\n(baseline)' = anno_points(tmp$cor,
                                                                                  pch = ifelse(tmp$p.adj < 0.05, 19, 1),
                                                                                  gp = gpar(col = ifelse(tmp$p.adj > 0.05, 'grey40',
                                                                                                         ifelse(tmp$cor<0, cols[1], cols[4])))))
    }

    
    p <- Heatmap(as.matrix(tmp[,c(2, 4, 6)]),
                 
                 row_labels = tmp$pop_name,
                 row_order = order(tmp$cor, decreasing = T),
                 row_title_rot = 0,
                 row_names_gp = grid::gpar(fontsize = 5.5),
                 row_title_gp = grid::gpar(fontsize = 10),
                 row_names_side = 'right',
                 row_split = tmp$type,
                 
                 cell_fun = function(j, i, x, y, width, height, fill) {
                   grid.text(as.matrix(tmp[,c(3,5,7)])[i, j], x, y, gp = gpar(fontsize = 10))
                 },
                 
                 cluster_columns = F,
                 column_title = NULL,
                 show_row_dend = F, show_column_dend = F,
                 
                 name = 'estimate\ncompliance\neffect\n(scaled)',
                 na_col = 'white',
                 
                 col = circlize::colorRamp2(breaks = seq(-1, 1, length.out = 30),
                                            colors = colorRampPalette(c(cols[c(1, 2)], "grey95", cols[c(5, 4)]))(30)
                                            #                            # colors = rev(color("batlow")(5)),
                                            #                            # colors = rev(cols[c(2,5,4, 6, 7)])
                 ),
                 column_names_gp = grid::gpar(fontsize = 9),
                 border_gp = gpar(lwd = 0.5),
                 border = T,
                 left_annotation = left_anno)  
    
  } else if(type == 'intervention'){
    
    tmp <- lmm_data |> 
      dplyr::select(name, 
                    `estimate_visitIdM2:interventionIdK`, `p.value_visitIdM2:interventionIdK`,
                    `estimate_visitIdM4:interventionIdK`, `p.value_visitIdM4:interventionIdK`,
                    `estimate_visitIdM6:interventionIdK`, `p.value_visitIdM6:interventionIdK`,
                    cor, p.adj,
                    pop_name, type) |> 
      dplyr::mutate(across(contains("p.value"), ~ ifelse(.<0.05, "*", ""))) |> 
      dplyr::rename_at(vars(contains('estimate_')), ~ gsub("estimate_visitId|:interventionIdK", "", .)) |> 
      dplyr::filter(!is.na(pop_name))
    
    if(!is.null(corr_baseline)){
      # annotation 
      left_anno <- rowAnnotation('Correlation with age\n(baseline)' = anno_points(tmp$cor,
                                                                                  pch = ifelse(tmp$p.adj < 0.05, 19, 1),
                                                                                  gp = gpar(col = ifelse(tmp$p.adj > 0.05, 'grey40',
                                                                                                         ifelse(tmp$cor<0, cols[1], cols[4])))))
    }
    
    p <- Heatmap(as.matrix(tmp[,c(2, 4, 6)]),
                 
                 row_labels = tmp$pop_name,
                 row_order = order(tmp$cor, decreasing = T),
                 row_title_rot = 0,
                 row_names_gp = grid::gpar(fontsize = 5.5),
                 row_title_gp = grid::gpar(fontsize = 10),
                 row_names_side = 'right',
                 row_split = tmp$type,
                 
                 cell_fun = function(j, i, x, y, width, height, fill) {
                   grid.text(as.matrix(tmp[,c(3,5,7)])[i, j], x, y, gp = gpar(fontsize = 10))
                 },
                 
                 cluster_columns = F,
                 column_title = NULL,
                 show_row_dend = F, show_column_dend = F,
                 
                 name = 'estimate\nMCT\neffect\n(scaled)',
                 na_col = 'white',
                 
                 col = circlize::colorRamp2(breaks = seq(-1, 1, length.out = 30),
                                            colors = colorRampPalette(c(cols[c(1, 2)], "grey95", cols[c(5, 4)]))(30)
                                            #                            # colors = rev(color("batlow")(5)),
                                            #                            # colors = rev(cols[c(2,5,4, 6, 7)])
                 ),
                 column_names_gp = grid::gpar(fontsize = 9),
                 border_gp = gpar(lwd = 0.5),
                 border = T,
                 left_annotation = left_anno)  
    
  }
  return(p)
  
}
