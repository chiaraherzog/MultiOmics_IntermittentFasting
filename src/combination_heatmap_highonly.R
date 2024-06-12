combination_heatmap_high <- function(lmm_data,
                                     cols = cols,
                                     relabel = labels,
                                     methyl_sub = NULL){
  
  tmp <- lmm_data |> 
    dplyr::select(x, 
                  `estimate_time`, `p.value_time`,
                  `estimate_time:interventionIdK`, `p.value_time:interventionIdK`) |> 
    dplyr::mutate(across(contains("p.value"), ~ ifelse(.<0.05, "*", ""))) |> 
    dplyr::rename(time = `estimate_time`,
                  `modulation;\nMCT` = `estimate_time:interventionIdK`)
  
  # condition: if assay name in x, need to remove this
  if(any(grepl("Composite", tmp$x))){
    
    tmp <- tmp |> 
      tidyr::separate(x, into = c('assay', "x"),
                      sep = "_", 
                      extra = 'merge')
    
    # methyl subset
    if(!is.null(methyl_sub)){
      tmp <- tmp |> 
        dplyr::filter(grepl(methyl_sub, assay) & grepl("methylation", assay))
      
    }
    
    tmp <- tmp |> 
      dplyr::select(-assay)
    
  }


  tmp <- tmp |> 
    # relabel
    dplyr::filter(x %in% relabel$x) |> 
    dplyr::left_join(relabel) |> 
    dplyr::mutate(x = label)
  
  ha = columnAnnotation(estimate = anno_mark(at = c(1, 2),
                                             labels = colnames(tmp[c(2, 4)]),
                                             which = 'column',
                                             labels_rot = 60,
                                             padding = 9,
                                             labels_gp = grid::gpar(fontsize = 10),
                                             link_height = unit(4, 'mm')))
  
  p <- Heatmap(as.matrix(tmp[,c(2, 4)]),
               row_labels = tmp$x,
               top_annotation = ha,
               cell_fun = function(j, i, x, y, width, height, fill) {
                 grid.text(as.matrix(tmp[,c(3, 5)])[i, j], x, y, gp = gpar(fontsize = 10))
               },
               row_names_side = 'left',
               column_names_side = 'top',
               # cluster_rows = T,
               row_order = order(tmp$time, decreasing = T),
               cluster_columns = F,
               column_title = NULL,
               show_row_dend = F, show_column_dend = F,
               show_column_names = F,
               name = 'estimate\n(scaled)',
               na_col = 'white',
               col = circlize::colorRamp2(breaks = seq(-max(abs(tmp[,c(2, 4)])), max(abs(tmp[,c(2, 4)])), length.out = 30),
                                          colors = colorRampPalette(c(cols[c(1, 2)], "grey95", cols[c(5, 4)]))(30)
                                          #                            # colors = rev(color("batlow")(5)),
                                          #                            # colors = rev(cols[c(2,5,4, 6, 7)])
               ),
               row_names_gp = grid::gpar(fontsize = 8),
               column_names_gp = grid::gpar(fontsize = 9),
               border_gp = gpar(lwd = 0.5),
               border = T) 
  
  return(p)
  
}

