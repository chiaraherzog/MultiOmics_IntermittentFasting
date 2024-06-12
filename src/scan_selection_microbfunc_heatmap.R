## scanning for enzymes in lmms run on functional predictions microbiomes
# plot heatmap in style of fragmented heatmap v2

scan_selection_microbfunc_heatmap <- function(exp,
                                              lmm_data_time,
                                              lmm_data_int,
                                              cols = c("#1b69a1",
                                                       "#48a0af",
                                                       "#f39668",
                                                       "#ec6669"),
                                              selection,
                                              pval=NULL){
  if(!require("gtools")){
    install.packages("gtools")
  }
  
  require(stringr)
  require(viridisLite)
  
  ## load output lmm models and keep only selection ##----
  
  tmp1 <- lmm_data_time %>%
    filter(str_detect(x, exp)) %>%
    dplyr::select(!contains(c("std.", "age_at_consent", "bmi_at_consent","interventionIdK")))
  
  tmp3 <- lmm_data_int %>%
    filter(str_detect(x, exp)) %>%
    dplyr::select(!contains(c("std.", "age_at_consent", "bmi_at_consent",
                              "estimate_interventionIdK",
                              "p.value_interventionIdK"))) %>% 
    dplyr::select(!ends_with(c("visitIdM2", "visitIdM4", "visitIdM6")))
  
  tmp <- tmp1 %>% 
    dplyr::full_join(tmp3, by = "x")
  
  tmp$x = gsub(paste0(exp,"_"),"",tmp$x)
  
  # filter out significant rows only, when requested
  
  if (!is.null(pval)){
    tmp <- tmp %>%
      filter(rowSums(across(starts_with("p.value"), ~ . < pval)) > 0)
  }
  
  # filter out only relevant rows, replace pval with * symbols and add info enzymes
  tmp <- tmp %>%
    dplyr::filter(x %in% c(selection$x)) %>% 
    dplyr::mutate(across(contains("p.value"), ~ ifelse(is.na(.), "", gtools::stars.pval(.)))) %>%
    dplyr::mutate(across(contains("p.value"), ~ gsub("[.]", "", .))) %>%
    dplyr::left_join(selection) %>%
    dplyr::mutate(x=label)
  
  # order
  tmp <- tmp %>% arrange(metabolite.class,Synthesis)
  tmp$metabolite.class <- factor(tmp$metabolite.class, levels = unique(tmp$metabolite.class))
  tmp$Synthesis <- factor(tmp$Synthesis, levels = unique(tmp$Synthesis))
  
  # Annotation
  cols_class = c("#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","#CC79A7","#FB28BE","#CF83E6","#92B4BA","#88DDD5")
  row_df <- tmp %>%
    select(label,Synthesis, metabolite.class) %>%
    column_to_rownames(var="label")
  row_ha = rowAnnotation(df = row_df,
                         col = list(Synthesis = 
                                      setNames(inferno(length(unique(tmp$Synthesis))),unique(tmp$Synthesis)),
                                    metabolite.class = 
                                      setNames(cols_class[1:length(unique(tmp$metabolite.class))],unique(tmp$metabolite.class))),
                         show_annotation_name = c(Synthesis = FALSE,metabolite.class = FALSE))
  
  ## Column splits and labels ##---- 
  col_split <- rep(c("time", "MCT"), each = 3)
  col_split <- factor(col_split, levels = c("time", "MCT"))
  col_labs <- rep(c("M2", "M4", "M6"), 2)
  
  # Indices
  tmp <- tmp %>% column_to_rownames(var = "x")
  ind_est <- grep("estimate", colnames(tmp))  
  ind_p <- grep("p.value", colnames(tmp))
  
  ## Draw heatmap ##---- 
  
  p <- Heatmap(as.matrix(tmp[,ind_est]),
               name = 'estimate\n(scaled)',
               
               # Row details
               row_names_side = 'left',
               #row_order = order(tmp$metabolite.class),
               cluster_rows = F,
               show_row_dend = F, 
               row_title = NULL,
               
               # Row annotation
               right_annotation = row_ha,
               
               # Column details
               column_labels = col_labs,
               column_title_gp = grid::gpar(fontsize = 10),
               column_split = col_split,
               cluster_columns = F,
               show_column_dend = F,
               
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
                                                                      "grey95", cols[c(3, 4)]))(30)),
               
               # Titles
               row_names_gp = grid::gpar(fontsize = 9),
               column_names_gp = grid::gpar(fontsize = 9),
               
               # Borders
               border_gp = gpar(lwd = 0.5),
               border = T)  
  
  return(p)
  
}
