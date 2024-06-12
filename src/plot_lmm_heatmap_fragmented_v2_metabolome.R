#' @param exp which experiment should be plotted
#' @param lmm_data_time Minimal model that provides time coefficient (no interaction) (out_lmm$`Minimal model`)
#' @param lmm_data_compliance Basic model that provides time*compliance interaction (out_lmm$`Basic model with BMI`)
#' @param lmm_data_int High compliance intervention output. (out_lmm$`Intervention (higher compliance)`)
#' @param cols Colour palette (leave blank for default colours [recommended], but can be altered if desired)
#' @param relabel (optional). Data frame supplying labels for the plot. x = variable, label = relabelled variable
#' @param age_cor should correlation with age be plotted on the side? default to F
#' @param bmi_cor should correlation with bmi be plotted on the side? default to F
#' @param FDR should p-values be corrected? default to F
#' @param cluster Should rows be clustered? If left blank, will set to 'default' and rows will be ordered by M6 time estimate. If 'cluster' (or other), rows will be clustered.
#' @return ComplexHeatmap object
#' 
# version made for metabolomes, also plotting compound class (if unwanted use microbiome function instead)
# option to plot age and/or bmi correlations at baseline
# option to FDR correct p-values, for each term (and category) separately 
# all terms with p.corrected < 0.05 are plotted.
# for uncorrected p-values, all rows with p < 0.01 are plotted

# main function
plot_lmm_heatmap_frag_v2_metab <- function(exp,
                                     lmm_data_time,
                                     lmm_data_int,
                                     cols = c("#1b69a1",
                                              "#48a0af",
                                              "#f39668",
                                              "#ec6669"),
                                     relabel = NULL,
                                     age_cor = F,
                                     bmi_cor = F,
                                     FDR = F,
                                     cluster = 'default'){
  
  if(!require("gtools")){
    install.packages("gtools")
  }
  
  require(stringr)
  
  source("src/FDR_correct.R")
  

  
  ## load output lmm models and correct p-values if requested, keep only significant results ##----
  
  if (FDR == F) {
    
    tmp1 <- lmm_data_time %>%
      filter(str_detect(x, exp)) %>%
      dplyr::select(!contains(c("std.", "age_at_consent", "bmi_at_consent")))
    
    tmp3 <- lmm_data_int %>%
      filter(str_detect(x, exp)) %>%
      dplyr::select(!contains(c("std.", "age_at_consent", "bmi_at_consent", "compliancemedium",
                                "estimate_interventionIdK",
                                "p.value_interventionIdK"))) %>% 
      dplyr::select(!ends_with(c("visitIdM2", "visitIdM4", "visitIdM6")))
    
    tmp <- tmp1 %>% 
      dplyr::full_join(tmp3, by = "x")
    
    tmp$x = gsub(paste0(exp,"_"),"",tmp$x)
    
    # filter out only significant rows and replace pval with * symbols
    tmp <- tmp %>%
      filter(rowSums(across(starts_with("p.value"), ~ . < 0.01), na.rm = TRUE) > 0) %>% 
      dplyr::mutate(across(contains("p.value"), ~ ifelse(is.na(.), "", gtools::stars.pval(.)))) |> 
      dplyr::mutate(across(contains("p.value"), ~ gsub("[.]", "", .)))
    
  } else if (FDR == T){
    
    # simple model FDR correct
    tmp1 <- lmm_data_time %>%
      filter(str_detect(x, exp))
    
    tmp1_p <- tmp1 %>%
      dplyr::select(x,contains("p.value")) %>%
      pivot_longer(cols = contains("p.value"),names_to = "estimate", values_to = "p_value") %>%
      FDR_correct() %>%
      pivot_wider(id_cols = x, names_from = estimate, values_from = p_value)
    
    tmp1 <- tmp1 %>%
      select(!contains("p.value")) %>%
      left_join(tmp1_p, by = "x") %>%
      dplyr::select(!contains(c("std.", "age_at_consent", "bmi_at_consent")))
    
    # MCT high compliance part FDR correct
    tmp3 <- lmm_data_int %>%
      filter(str_detect(x, exp))
    
    tmp3_p <- tmp3 %>%
      dplyr::select(x,contains("p.value")) %>%
      pivot_longer(cols = contains("p.value"),names_to = "estimate", values_to = "p_value") %>%
      FDR_correct() %>%
      pivot_wider(id_cols = x, names_from = estimate, values_from = p_value)
    
    tmp3 <- tmp3 %>%
      select(!contains("p.value")) %>%
      left_join(tmp3_p, by = "x") %>%
      dplyr::select(!contains(c("std.", "age_at_consent", "bmi_at_consent", "compliancemedium",
                                "estimate_interventionIdK",
                                "p.value_interventionIdK"))) %>% 
      dplyr::select(!ends_with(c("visitIdM2", "visitIdM4", "visitIdM6")))
    
    # combine
    tmp <- tmp1 %>% 
      dplyr::full_join(tmp3, by = "x")
    
    tmp$x = gsub(paste0(exp,"_"),"",tmp$x)
    
    # filter out only significant rows and replace pval with * symbols
    tmp <- tmp %>%
      filter(rowSums(across(starts_with("p.value"), ~ . < 0.05)) > 0, na.rm = TRUE) %>% 
      dplyr::mutate(across(contains("p.value"), ~ ifelse(is.na(.), "", gtools::stars.pval(.)))) |> 
      dplyr::mutate(across(contains("p.value"), ~ gsub("[.]", "", .)))
  }
  

  
  ## Class annotation ##----
  # Load the correct reference map  

  if (grepl("Urine", exp)) {
    RefMet <-  readRDS("out/RefMet_mapped_urine.Rds")
  } else if (grepl("Saliva", exp)) {
    RefMet <-  readRDS("out/RefMet_mapped_saliva.Rds")
  }
  
  # Harmonize input names
  tmp$x <- sub(paste0("^",exp,"_"), "", tmp$x)
  
  tmp$x <- gsub("^X\\d+\\.", "", tmp$x)
  
  # Metabolites that are unknown to RefMet - try more standardized name
  tmp$x[tmp$x=="Acetate.mM."] <- "Acetate"
  tmp$x[tmp$x=="Methyl.2.oxovalerate"] <- "3-Methyl-2-oxovaleric acid"
  tmp$x[tmp$x=="Aminopentanoate"] <- "5-Aminopentanoate"
  tmp$x[tmp$x=="TMA..N.oxide"] <- "Trimethylamine N-oxide"
  tmp$x[tmp$x=="X2.Hydroxyisobutyrate"] <- "2-hydroxybutyric acid"
  tmp$x[tmp$x=="X3.Hydroxyisobutyrate"] <- "3-hydroxybutyric acid"
  

  
  # Replace all "-" entries with "unknown" in merged_data dataframe
  # RefMet <- RefMet %>%
  #   mutate_all(~replace(., . == NA, "unknown"))
  
  
  tmp <- merge(RefMet, tmp, by.x = "Input.name", by.y = "x", all.x = F, all.y = TRUE)
  tmp$Standardized.name <- NULL
  tmp$Formula <- NULL
  tmp$x <- tmp$Input.name

  
  cols_class <- c("Nicotinic acid alkaloids" = "#1b69a1", "Amino acids and peptides" = "#48a0af", "-" = "#71b5a9", 
            "Tryptophan alkaloids" = "#ec6669", "Phenolic acids" = "#f39668", "Fatty acids" = "#bd647d", 
            "Carbonyl compounds" = "#395262", "Monosaccharides" = "#5f70a8", "Cholines" = "#d4a34b", "TCA acids" = "#4da457", 
            "Azoles" = "#d16fa8", "Organonitrogen compounds" = "#78b646", "Carboxylic acids" = "#3b738f", 
            "Short-chain acids" = "#a36f2d", "Sulfonic acids" = "#d4a34b", "Amines" = "#6e2b57", "Amine oxides" = "#c97f2e", 
            "Organic carbonic acids" = "#8e4555", "Hydrocarbons" = "#4f6d32", "Phenylpropanoids" = "#ec6669", 
            "Primary alcohols" = "#9c572c", "Alcohols and polyols" = "#5a9d79", "Phenols" = "#d09636", 
            "Fatty amines" = "#7d3b65", "Pyrimidines" = "#1b69a1", "Purines" = "#48a0af", "Keto acids"="#832c9b")
  
  
  row_ha = rowAnnotation(Main.class = tmp$Main.class,
                         col = list(Main.class = cols_class),
                         show_annotation_name = c(Main.class = FALSE))
  
  ## Column splits and labels ##---- 
  col_split <- rep(c("time", "MCT"), each = 3)
  col_split <- factor(col_split, levels = c("time", "MCT"))
  col_labs <- rep(c("M2", "M4", "M6"), 2)
  
  if(!is.null(relabel)){
    
    tmp <- tmp |> 
      # relabel
      dplyr::filter(x %in% relabel$x) |> 
      dplyr::left_join(relabel) |> 
      dplyr::mutate(x = label)
    
  }
  
  row_labs <- ifelse(tmp$`p.value_visitIdM6:interventionIdK` == "*" |
                       tmp$`p.value_visitIdM6:interventionIdK` == "**" |
                       tmp$`p.value_visitIdM6:interventionIdK` == "***" |
                       tmp$`p.value_visitIdM6:interventionIdK` == ".",
                     paste0("**", tmp$x, "**"),
                     tmp$x)
  
  # Indices
  ind_est <- grep("estimate", colnames(tmp))  
  ind_p <- grep("p.value", colnames(tmp))

  ## Draw heatmap ##---- 
  
  if(cluster == 'default'){
    p <- Heatmap(as.matrix(tmp[,ind_est]),
                 name = 'estimate\n(scaled)',
                 
                 # Row details
                 row_labels = gt_render(row_labs),
                 row_names_side = 'left',
                 row_order = order(tmp$estimate_visitIdM6, decreasing = T),
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
  } else {
    p <- Heatmap(as.matrix(tmp[,ind_est]),
                 #name = 'estimate\n(scaled)',
                 
                 # Row details
                 row_labels = gt_render(row_labs),
                 row_names_side = 'left',
                 cluster_rows = T,
                 cluster_row_slices = F,
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
  }
  
  return(p)
  
}
