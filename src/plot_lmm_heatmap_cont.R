plot_lmm_heatmap_cont <- function(lmm_data,
                             cols = cols,
                             type = 'compliance'){
  
  if(type == 'compliance'){
    
    tmp <- lmm_data |> 
      dplyr::select(x, 
                    `estimate_time:compliancehigh`, `p.value_time:compliancehigh`) |> 
      dplyr::mutate(across(contains("p.value"), ~ ifelse(.<0.05, "*", ""))) |> 
      dplyr::rename_at(vars(contains('estimate_')), ~ gsub("estimate_time:", "", .))
    
    p <- Heatmap(as.matrix(tmp[,2]),
                 row_labels = tmp$x,
                 cell_fun = function(j, i, x, y, width, height, fill) {
                   grid.text(as.matrix(tmp[,c(3)])[i, j], x, y, gp = gpar(fontsize = 10))
                 },
                 row_names_side = 'left',
                 cluster_rows = T,
                 cluster_columns = F,
                 column_title = NULL,
                 show_row_dend = F, show_column_dend = F,
                 name = 'estimate\ncompliance\neffect\n(scaled)',
                 na_col = 'white',
                 col = circlize::colorRamp2(breaks = seq(-max(abs(tmp[,2])), max(abs(tmp[,2])), length.out = 30),
                                            colors = colorRampPalette(c(cols[c(1, 2)], "grey95", cols[c(5, 4)]))(30)
                                            #                            # colors = rev(color("batlow")(5)),
                                            #                            # colors = rev(cols[c(2,5,4, 6, 7)])
                 ),
                 row_names_gp = grid::gpar(fontsize = 9),
                 column_names_gp = grid::gpar(fontsize = 9),
                 border_gp = gpar(lwd = 0.5),
                 border = T)  
    
  } else if(type == 'intervention'){
    
    tmp <- lmm_data |> 
      dplyr::select(x, 
                    `estimate_time:interventionIdK`, `p.value_time:interventionIdK`) |> 
      dplyr::mutate(across(contains("p.value"), ~ ifelse(.<0.05, "*", ""))) |> 
      dplyr::rename_at(vars(contains('estimate_')), ~ gsub("estimate_time:|:interventionIdK", "", .))
    
    p <- Heatmap(as.matrix(tmp[,c(2)]),
                 row_labels = tmp$x,
                 cell_fun = function(j, i, x, y, width, height, fill) {
                   grid.text(as.matrix(tmp[,c(3)])[i, j], x, y, gp = gpar(fontsize = 10))
                 },
                 row_names_side = 'left',
                 cluster_rows = T,
                 cluster_columns = F,
                 column_title = NULL,
                 show_row_dend = F, show_column_dend = F,
                 name = 'estimate\nintervention\neffect\n(scaled)',
                 na_col = 'white',
                 col = circlize::colorRamp2(breaks = seq(-max(abs(tmp[,2])), max(abs(tmp[,2])), length.out = 30),
                                            colors = colorRampPalette(c(cols[c(1, 2)], "grey95", cols[c(5, 4)]))(30)
                                            #                            # colors = rev(color("batlow")(5)),
                                            #                            # colors = rev(cols[c(2,5,4, 6, 7)])
                 ),
                 row_names_gp = grid::gpar(fontsize = 9),
                 column_names_gp = grid::gpar(fontsize = 9),
                 border_gp = gpar(lwd = 0.5),
                 border = T)  
    
  }
  return(p)
  
}
