combination_heatmap_high_factor <- function(lmm_data,
                                     cols = cols,
                                     relabel = NULL,
                                     methyl_sub = NULL,
                                     m6_only = F,
                                     m6_sep = F){
  
  tmp <- lmm_data |> 
    dplyr::select(x, 
                  `estimate_visitIdM2`, `p.value_visitIdM2`,
                  `estimate_visitIdM2:interventionIdK`, `p.value_visitIdM2:interventionIdK`,
                  `estimate_visitIdM4`, `p.value_visitIdM4`,
                  `estimate_visitIdM4:interventionIdK`, `p.value_visitIdM4:interventionIdK`,
                  `estimate_visitIdM6`, `p.value_visitIdM6`,
                  `estimate_visitIdM6:interventionIdK`, `p.value_visitIdM6:interventionIdK`) |> 
    dplyr::mutate(across(contains("p.value"), ~ ifelse(.<0.05 & !is.na(.), "*", "")))
  
  # condition: if assay name in x, need to remove this
  if(any(grepl("Composite|Functional|Blood haemogram", tmp$x))){
    
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
  
  if(!is.null(relabel)){
    
  tmp <- tmp |> 
    # relabel
    dplyr::filter(x %in% relabel$x) |> 
    dplyr::left_join(relabel) |> 
    dplyr::mutate(x = label)
  
  }
  

  
  if(m6_only == F){
    
    # anyMissing
    
    col_split = c(rep("M2", 2), rep("M4", 2), rep("M6", 2))
    
    if(m6_sep == T){
    row_split = ifelse(is.na(tmp$`estimate_visitIdM2:interventionIdK`), 2, 1)
    } else {
      row_split = rep(1, nrow(tmp))
    }
    
    ha = columnAnnotation(estimate = rep(c("time", "modulation MCT"), 3),
                          col = list(estimate = c('time' = cols[3],
                                                  'modulation MCT' = cols[8])))
    
  p <- Heatmap(as.matrix(tmp[,c(2, 4, 6, 8, 10, 12)]),
               row_labels = tmp$x,
               top_annotation = ha,
               column_split = col_split,
               row_split = row_split,
               cell_fun = function(j, i, x, y, width, height, fill) {
                 grid.text(as.matrix(tmp[,c(3, 5, 7, 9, 11, 13)])[i, j], x, y, gp = gpar(fontsize = 10))
               },
               row_names_side = 'left',
               column_names_side = 'top',
               # cluster_rows = T,
               row_order = order(tmp$estimate_visitIdM6, decreasing = T),
               cluster_columns = F,
               # column_title = NULL,
               show_row_dend = F, show_column_dend = F,
               show_column_names = F,
               row_title = NULL,
               name = 'estimate\n(scaled)',
               na_col = 'white',
               col = circlize::colorRamp2(breaks = seq(-max(abs(na.omit(tmp[,c(2, 4,6,8,10,12)]))), max(abs(na.omit(tmp[,c(2, 4,6,8,10,12)]))), length.out = 30),
                                          colors = colorRampPalette(c(cols[c(1, 2)], "grey95", cols[c(5, 4)]))(30)
                                          #                            # colors = rev(color("batlow")(5)),
                                          #                            # colors = rev(cols[c(2,5,4, 6, 7)])
               ),
               row_names_gp = grid::gpar(fontsize = 9),
               column_names_gp = grid::gpar(fontsize = 9),
               border_gp = gpar(lwd = 0.5),
               border = T) 
  } else {
    
    col_split = c(rep("M6", 2))
    
    ha = columnAnnotation(estimate = c("time", "modulation MCT"),
                          col = list(estimate = c('time' = cols[3],
                                                  'modulation MCT' = cols[8])))
    
    p <- Heatmap(as.matrix(tmp[,c(10, 12)]),
                 row_labels = tmp$x,
                 top_annotation = ha,
                 column_split = col_split,
                 cell_fun = function(j, i, x, y, width, height, fill) {
                   grid.text(as.matrix(tmp[,c(11, 13)])[i, j], x, y, gp = gpar(fontsize = 10))
                 },
                 row_names_side = 'left',
                 column_names_side = 'top',
                 # cluster_rows = T,
                 row_order = order(tmp$estimate_visitIdM6, decreasing = T),
                 cluster_columns = F,
                 # column_title = NULL,
                 show_row_dend = F, show_column_dend = F,
                 show_column_names = F,
                 name = 'estimate\n(scaled)',
                 na_col = 'white',
                 col = circlize::colorRamp2(breaks = seq(-max(abs(na.omit(tmp[,c(10,12)]))), max(abs(na.omit(tmp[,c(10,12)]))), length.out = 30),
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

