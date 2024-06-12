signASV_lmm_custom_comb <- function(exp,
                                    sample_type,
                                    lmm_data,
                                    estimates,
                                    datASV,
                                    label_types,
                                    FDR = T,
                                    mean = T,
                                    visitId = NULL){
  
  source("src/FDR_correct.R")
  
  dat2 <- list()
  
  for (i in 1:length(exp)){
  
    if (FDR == F) {
   
        tmp = lmm_data |> 
          filter(str_detect(x, exp[i])) %>%
          dplyr::select(x, paste0("p.value_",estimates)) 
        
        # rename features
        tmp$x = gsub(paste0(exp[i],"_"),"",tmp$x)
        
        # filter out only significant rows
        tmp = tmp %>%
          filter(rowSums(across(starts_with("p.value") & (contains("visitId")|contains("interventionId")), ~ . < 0.01)) > 0) %>% 
          dplyr::mutate(across(contains("p.value"), ~ ifelse(is.na(.), "", gtools::stars.pval(.)))) |> 
          dplyr::mutate(across(contains("p.value"), ~ gsub("[.]", "", .)))
 
    } else if (FDR == T){
      
      # FDR correction
      tmp = lmm_data |> 
        filter(str_detect(x, exp[i])) %>%
        dplyr::select(x,contains("p.value")) %>%
        pivot_longer(cols = contains("p.value"),names_to = "estimate", values_to = "p_value") %>%
        FDR_correct() %>%
        pivot_wider(id_cols = x, names_from = estimate, values_from = p_value) %>%
        dplyr::select(x,paste0("p.value_",estimates)) 
      
      # rename features
      tmp$x = gsub(paste0(exp[i],"_"),"",tmp$x)
      
      # filter out only significant rows
      tmp <- tmp %>%
        filter(rowSums(across(starts_with("p.value") & (contains("visitId")|contains("interventionId")), ~ . < 0.05)) > 0) %>% 
        dplyr::mutate(across(contains("p.value"), ~ ifelse(is.na(.), "", gtools::stars.pval(.)))) |> 
        dplyr::mutate(across(contains("p.value"), ~ gsub("[.]", "", .)))
    }
  
    # Relative abundance
    dat = datASV[[i]] %>% 
      filter(ASV %in% tmp$x) %>%
      mutate(label = paste(Family,ASV, sep = ";")) %>%
      select(label,contains(sample_type[i]))
    
    if (mean == T){
      
      # calculate mean
      dat = dat %>%
        pivot_longer(cols = contains(sample_type[i]), names_to = "sampleId", values_to = "rel_change_abund") %>%
        mutate(subjectId = substr(sampleId, 1, 4)) %>%
        mutate(visitId = substr(sampleId, 5, 6)) %>%
        select(-sampleId) %>%
        group_by(label, visitId) %>%
        summarize(mean_rel_change_abund = mean(rel_change_abund, na.rm=TRUE)) %>%
        pivot_wider(names_from = visitId, values_from = mean_rel_change_abund)
      
      # add p-values and phylum info
      dat = dat %>%
        mutate(ASV = gsub(".*?;", "", label)) %>%
        left_join(tmp, by = c("ASV"="x")) %>%
        left_join(datASV[[i]][c("ASV", "Phylum")]) 
      
      dat$sample_type <- label_types[i]
      dat2[[i]] <- dat
      
    } else if (mean ==F){
      
      # add phylum info
      dat = dat %>%
        mutate(ASV = gsub(".*?;", "", label)) %>%
        left_join(datASV[[i]][c("ASV", "Phylum")]) %>%
        select(label,ASV,Phylum,contains(visitId)) # plot only for one timepoint at a time
      
      dat <- dat[, order(names(dat))]
      
      # keep only participant + visitId in colnames so can be combined, add sample type designation
      colnames(dat) <- gsub(sample_type[i],"",colnames(dat))
      dat$sample_type <- label_types[i]

      dat2[[i]] <- dat
  
    }
  }
  
  # combine dat list and make final heatmap
  
  dat2 <- do.call(full_join, dat2) %>%
    column_to_rownames(var="label")
  
  # final heatmap
  
  if (mean == T){
    
    # final heatmap
    ind_main <- which(colnames(dat) %in% c("M2","M4","M6")) 
    ind_p <- grep("p.value", colnames(dat))
    
    # add row_ha ~ sample type
    
    cols = c("#1b69a1","#48a0af", "#f39668","#ec6669")
    p <- Heatmap(as.matrix(dat[,ind_main]),
                 name = 'mean change from baseline\n(relative abundance)',
                 
                 # Row details
                 row_names_side = 'left',
                 row_order = order(dat[,ind_main]$M6, decreasing = T),
                 show_row_dend = F, 
                 row_title = NULL,
                 
                 # Column details
                 cluster_columns = F,
                 show_column_dend = F,
                 
                 # Annotation
                 cell_fun = function(j, i, x, y, width, height, fill) {
                   grid.text(as.matrix(dat[,ind_p])[i, j], x, y, gp = gpar(fontsize = 7), just = "centre")
                 },
                 
                 # Colours
                 na_col = 'white',
                 col = circlize::colorRamp2(breaks = seq(-max(abs(na.omit(dat[,ind_main]))),
                                                         max(abs(na.omit(dat[,ind_main]))),
                                                         length.out = 30),
                                            colors = colorRampPalette(c(cols[c(1, 2)],
                                                                        "grey95", cols[c(3, 4)]))(30)
                 ),
                 
                 # Titles
                 row_names_gp = grid::gpar(fontsize = 9),
                 column_names_gp = grid::gpar(fontsize = 9),
                 
                 # Borders
                 border_gp = gpar(lwd = 0.5),
                 border = T) 
    
  } else if (mean ==F){
    
    ind_main <- grep(visitId, colnames(dat2))
    
    # column and row annotation
    column_groups <- substr(colnames(dat2[,ind_main]), 1, 1)
    column_ha = HeatmapAnnotation(InterventionId = column_groups,
                                  col = list(InterventionId = c("I" = "#5DBAB2", "K" = "#832c9b")),
                                  show_annotation_name = c(InterventionId = FALSE))
    
    # add row_ha
    cols_type = c("#48a0af","#964B00") #set for saliva, stool in microbiome
    cols_type = setNames(cols_type, label_types)
    # dat$Phylum_Colors = cols_phyla[dat$Phylum] # do I need this??
    dat2$sample_type = factor(dat2$sample_type, levels = label_types)
    row_ha = rowAnnotation(sample_type = dat2$sample_type,
                           col = list(sample_type = cols_type),
                           show_annotation_name = c(sample_type = FALSE))
    
    cols = c("#1b69a1","#48a0af", "#f39668","#ec6669")
    p <- Heatmap(as.matrix(dat2[,ind_main]),
                 name = paste0(visitId, '\nchange from baseline\n(rel. abund.)'),
                 
                 # Row details
                 row_names_side = 'left',
                 cluster_rows = T,
                 show_row_dend = F, 
                 row_title = NULL,
                 
                 # Column details
                 cluster_columns = F,
                 show_column_dend = F,
                 
                 # Annotation
                 top_annotation = column_ha,
                 right_annotation = row_ha,
                 
                 # Colours
                 na_col = 'white',
                 col = circlize::colorRamp2(breaks = seq(-max(abs(dat2[,ind_main]), na.rm = TRUE),
                                                         max(abs(dat2[,ind_main]), na.rm = TRUE),
                                                         length.out = 30),
                                            colors = colorRampPalette(c(cols[c(1, 2)],
                                                                        "grey95", cols[c(3, 4)]))(30)
                 ),
                 
                 # Titles
                 row_names_gp = grid::gpar(fontsize = 9),
                 show_column_names = FALSE,
                 
                 # Borders
                 border_gp = gpar(lwd = 0.5),
                 border = T) 
    
  }
  
  
  return(p)
}