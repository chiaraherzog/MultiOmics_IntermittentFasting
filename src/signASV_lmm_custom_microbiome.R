signASV_lmm_custom <- function(exp,
                               sample_type,
                               lmm_data,
                               estimates,
                               datASV,
                               buccal_ic_cor = F,
                               FDR = T,
                               mean = T,
                               visitId = NULL){
  
  source("src/FDR_correct.R")
  
  if (FDR == F) {
    
    tmp = lmm_data |> 
      filter(str_detect(x, exp)) %>%
      dplyr::select(x, paste0("p.value_",estimates)) 
    
    # rename features
    tmp$x = gsub(paste0(exp,"_"),"",tmp$x)
    
    # filter out only significant rows
    tmp = tmp %>%
      filter(rowSums(across(starts_with("p.value") & (contains("visitId")|contains("interventionId")), ~ . < 0.01)) > 0) %>% 
      dplyr::mutate(across(contains("p.value"), ~ ifelse(is.na(.), "", gtools::stars.pval(.)))) |> 
      dplyr::mutate(across(contains("p.value"), ~ gsub("[.]", "", .)))
    
  } else if (FDR == T){
    
    # FDR correction
    tmp = lmm_data |> 
      filter(str_detect(x, exp)) %>%
      dplyr::select(x,contains("p.value")) %>%
      pivot_longer(cols = contains("p.value"),names_to = "estimate", values_to = "p_value") %>%
      FDR_correct() %>%
      pivot_wider(id_cols = x, names_from = estimate, values_from = p_value) %>%
      dplyr::select(x,paste0("p.value_",estimates)) 
    
    # rename features
    tmp$x = gsub(paste0(exp,"_"),"",tmp$x)
    
    # filter out only significant rows
    tmp <- tmp %>%
      filter(rowSums(across(starts_with("p.value") & (contains("visitId")|contains("interventionId")), ~ . < 0.05)) > 0) %>% 
      dplyr::mutate(across(contains("p.value"), ~ ifelse(is.na(.), "", gtools::stars.pval(.)))) |> 
      dplyr::mutate(across(contains("p.value"), ~ gsub("[.]", "", .)))
    
  }
  
  # Relative abundance
  dat = datASV %>% 
    filter(ASV %in% tmp$x) %>%
    mutate(label = paste(Family,ASV, sep = ";")) %>%
    select(label,contains(sample_type))
  
  if (mean == T){
    
    # calculate mean
    dat = dat %>%
      pivot_longer(cols = contains(sample_type), names_to = "sampleId", values_to = "rel_change_abund") %>%
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
      left_join(datASV[c("ASV", "Phylum")]) 
    
    # Baseline correlation
    if(buccal_ic_cor == T) {
      
      load("out/corrIcBaseline_saliva_raw.R")
      corr <- corr |> 
        dplyr::rename(ASV=rowname) |> 
        dplyr::slice(match(dat$ASV, ASV))
      
      dat <- dat |> 
        dplyr::left_join(dplyr::select(corr, ASV, cor, p), by = c('ASV' = 'ASV')) %>%
        column_to_rownames(var="label")
      
      hr <- rowAnnotation('buccal_ic\n(baseline)' = anno_points(dat$cor,
                                                                pch = ifelse(dat$p < 0.05, 19, 1),
                                                                ylim = c(-1, 1),
                                                                gp = gpar(col = ifelse(dat$cor<0, cols[1], cols[4]))))
    } else{
      dat <- dat %>%
        column_to_rownames(var="label")
      hr <- NULL
    }
    
    # final heatmap
    ind_main <- which(colnames(dat) %in% c("M2","M4","M6")) 
    ind_p <- grep("p.value", colnames(dat))
    
    cols = c("#1b69a1","#48a0af", "#f39668","#ec6669")
    p <- Heatmap(as.matrix(dat[,ind_main]),
                 name = 'mean change from baseline\n(relative abundance)',
                 
                 # Row details
                 row_names_side = 'left',
                 row_order = order(dat[,ind_main]$M6, decreasing = T),
                 show_row_dend = F, 
                 row_title = NULL,
                 
                 # Row annotation
                 right_annotation = hr,
                 
                 # Column details
                 cluster_columns = F,
                 show_column_dend = F,
                 
                 # Annotation
                 cell_fun = function(j, i, x, y, width, height, fill) {
                   grid.text(as.matrix(dat[,ind_p])[i, j], x, y, gp = gpar(fontsize = 7), just = "centre")
                 },
                 
                 # Colours
                 na_col = 'white',
                 col = circlize::colorRamp2(breaks = seq(-max(abs(dat[,ind_main]), na.rm = TRUE),
                                                         max(abs(dat[,ind_main]), na.rm = TRUE),
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
    
    # add phylum info
    dat = dat %>%
      mutate(ASV = gsub(".*?;", "", label)) %>%
      left_join(datASV[c("ASV", "Phylum")]) %>%
      #column_to_rownames(var="label") %>%
      select(label,ASV,Phylum,contains(visitId)) # plot only for one timepoint at a time

    # row annotation, including baseline correlation with ic when requested
    if(buccal_ic_cor == T) {
      
      load("out/corrIcBaseline_saliva_raw.R")
      corr <- corr |> 
        dplyr::rename(ASV=rowname) |> 
        dplyr::slice(match(dat$ASV, ASV))
      
      dat <- dat |> 
        dplyr::left_join(dplyr::select(corr, ASV, cor, p), by = c('ASV' = 'ASV')) %>%
        column_to_rownames(var="label")
      
      cols_phyla = c("#D81B60","#1E88E5","#FFC107","#004D40","#53A20A","#B35F50","#97027B","#559354","#5434CA","#BFFB22")
      cols_phyla = cols_phyla[1:length(unique(dat$Phylum))]
      cols_phyla = setNames(cols_phyla, unique(dat$Phylum))
      dat$Phylum_Colors = cols_phyla[dat$Phylum] # do I need this??
      dat$Phylum = factor(dat$Phylum, levels = unique(dat$Phylum))
      row_ha = rowAnnotation(df = dat$Phylum,
                             col = list(Phylum = cols_phyla),
                             #show_annotation_name = c(Phylum = FALSE),
                             'buccal_ic\n(baseline)' = anno_points(dat$cor,
                                                                   pch = ifelse(dat$p < 0.05, 19, 1),
                                                                   ylim = c(-1, 1),
                                                                   gp = gpar(col = ifelse(dat$cor<0, cols[1], cols[4]))))
      
      dat <- dat[, order(names(dat))]
      
      ind_main <- grep(sample_type, colnames(dat))
      
      # column annotation
      column_groups <- substr(colnames(dat[,ind_main]), 1, 1)
      column_ha = HeatmapAnnotation(InterventionId = column_groups,
                                    col = list(InterventionId = c("I" = "#5DBAB2", "K" = "#832c9b")),
                                    show_annotation_name = c(InterventionId = FALSE))
      
      
    } else{
      dat = dat %>%
        column_to_rownames(var="label")
      
      dat <- dat[, order(names(dat))]
      
      ind_main <- grep(sample_type, colnames(dat))
      
      # column annotation
      column_groups <- substr(colnames(dat[,ind_main]), 1, 1)
      column_ha = HeatmapAnnotation(InterventionId = column_groups,
                                    col = list(InterventionId = c("I" = "#5DBAB2", "K" = "#832c9b")),
                                    show_annotation_name = c(InterventionId = FALSE))
      
      cols_phyla = c("#D81B60","#1E88E5","#FFC107","#004D40","#53A20A","#B35F50","#97027B","#559354","#5434CA","#BFFB22")
      cols_phyla = cols_phyla[1:length(unique(dat$Phylum))]
      cols_phyla = setNames(cols_phyla, unique(dat$Phylum))
      dat$Phylum_Colors = cols_phyla[dat$Phylum] # do I need this??
      dat$Phylum = factor(dat$Phylum, levels = unique(dat$Phylum))
      row_ha = rowAnnotation(Phylum = dat$Phylum,
                             col = list(Phylum = cols_phyla),
                             show_annotation_name = c(Phylum = FALSE))
    }

    # final heatmap
    #ind_p <- which(colnames(dat) %in% c("Phylum")) 
    
    cols = c("#1b69a1","#48a0af", "#f39668","#ec6669")
    p <- Heatmap(as.matrix(dat[,ind_main]),
                 name = paste0(visitId, ' change from baseline\n(relative abundance)'),
                 
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
                 col = circlize::colorRamp2(breaks = seq(-max(abs(dat[,ind_main]), na.rm = TRUE),
                                                         max(abs(dat[,ind_main]), na.rm = TRUE),
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