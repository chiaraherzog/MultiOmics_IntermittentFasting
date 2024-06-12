plot_lmm_heatmap <- function(lmm_data,
                             cols = cols,
                             type = 'compliance',
                             relabel = NULL,
                             filter_relabel = T,
                             methyl_sub = NULL){
  
  if(type == 'compliance'){
    
  tmp <- lmm_data |> 
    dplyr::select(x, 
                  `estimate_visitIdM2:compliancehigh`, `p.value_visitIdM2:compliancehigh`,
                  `estimate_visitIdM4:compliancehigh`, `p.value_visitIdM4:compliancehigh`,
                  `estimate_visitIdM6:compliancehigh`, `p.value_visitIdM6:compliancehigh`) |> 
    dplyr::mutate(across(contains("p.value"), ~ ifelse(.<0.05, "*", ""))) |> 
    dplyr::rename_at(vars(contains('estimate_')), ~ gsub("estimate_visitId|:compliancehigh", "", .))
  
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
  
  if(!is.null(relabel)){
    
    if(filter_relabel == T){
      tmp <- tmp |> 
      dplyr::filter(x %in% relabel$x)
    }
    
    tmp <- tmp |> 
      dplyr::left_join(relabel) |> 
      dplyr::mutate(x = label) |> 
      dplyr::select(-label)
  }

  p <- Heatmap(as.matrix(tmp[,c(2, 4, 6)]),
          row_labels = tmp$x,
          cell_fun = function(j, i, x, y, width, height, fill) {
            grid.text(as.matrix(tmp[,c(3,5,7)])[i, j], x, y, gp = gpar(fontsize = 10))
          },
          row_names_side = 'left',
          cluster_rows = T,
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
          row_names_gp = grid::gpar(fontsize = 9),
          column_names_gp = grid::gpar(fontsize = 9),
          border_gp = gpar(lwd = 0.5),
          border = T)  
  
  } else if(type == 'intervention'){
    
    tmp <- lmm_data |> 
      dplyr::select(x, 
                    `estimate_visitIdM2:interventionIdK`, `p.value_visitIdM2:interventionIdK`,
                    `estimate_visitIdM4:interventionIdK`, `p.value_visitIdM4:interventionIdK`,
                    `estimate_visitIdM6:interventionIdK`, `p.value_visitIdM6:interventionIdK`) |> 
      dplyr::mutate(across(contains("p.value"), ~ ifelse(.<0.05, "*", ""))) |> 
      dplyr::rename_at(vars(contains('estimate_')), ~ gsub("estimate_visitId|:interventionIdK", "", .))
    
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
    
    if(!is.null(relabel)){
      
      if(filter_relabel == T){
        tmp <- tmp |> 
          dplyr::filter(x %in% relabel$x)
      }
      
      tmp <- tmp |> 
        dplyr::left_join(relabel) |> 
        dplyr::mutate(x = label) |> 
        dplyr::select(-label)
    }
    
    p <- Heatmap(as.matrix(tmp[,c(2, 4, 6)]),
                 row_labels = tmp$x,
                 cell_fun = function(j, i, x, y, width, height, fill) {
                   grid.text(as.matrix(tmp[,c(3,5,7)])[i, j], x, y, gp = gpar(fontsize = 10))
                 },
                 row_names_side = 'left',
                 cluster_rows = T,
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
                 row_names_gp = grid::gpar(fontsize = 9),
                 column_names_gp = grid::gpar(fontsize = 9),
                 border_gp = gpar(lwd = 0.5),
                 border = T)  
    
  } 
    
  return(p)
  
}
