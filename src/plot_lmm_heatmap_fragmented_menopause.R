#' @param exp which experiment should be plotted
#' @param lmm_data_overall General model with menopauze term (no interaction) (out_lmm$`Menopause general`)
#' @param lmm_data_int Data from model with additional interaction term. (out_lmm$`Menopause interaction`)
#' @param cols Colour palette (leave blank for default colours [recommended], but can be altered if desired)
#' @param FDR should p-values be corrected? default to F
#' @return ComplexHeatmap object
#' 
# version made for additional checks on menopauzal status
# models run for high compliance individuals only
# option to FDR correct p-values, for each term (and category) separately 
# all terms with p.corrected < 0.05 are plotted.
# for uncorrected p-values, all rows with p < 0.05 are plotted

# main function
plot_lmm_heatmap_frag <- function(exp,
                                  lmm_data_overall,
                                  lmm_data_int,
                                  cols = c("#1b69a1",
                                           "#48a0af",
                                           "#f39668",
                                           "#ec6669"),
                                  FDR = F){
  
  if(!require("gtools")){
    install.packages("gtools")
  }
  
  library(stringr)
  library(purrr)
  library(openxlsx)
  library(rcartocolor)
  
  source("src/FDR_correct.R")
  
  ## load output lmm models and correct p-values if requested, keep only significant results ##----
  
  if (FDR == F) {
    
    tmp1 <- lmm_data_overall %>%
      dplyr::filter(purrr::map_lgl(x, ~ any(stringr::str_detect(.x, exp)))) %>%
      dplyr::select(x|contains(c("mpstatrsyes"))) %>%
      dplyr::select(!contains(c("std.")))
    
    tmp3 <- lmm_data_int %>%
      dplyr::filter(purrr::map_lgl(x, ~ any(stringr::str_detect(.x, exp)))) %>%
      dplyr::select(!contains(c("std.", "age_at_consent", "bmi_at_consent","interventionIdK",
                                "estimate_mpstatrsyes",
                                "p.value_mpstatrsyes"))) %>% 
      dplyr::select(!ends_with(c("visitIdM2", "visitIdM4", "visitIdM6")))
    
    tmp <- tmp1 %>% 
      dplyr::full_join(tmp3, by = "x")
    
    #tmp$x = gsub(paste0(exp,"_"),"",tmp$x)
    
    # filter out only significant rows and replace pval with * symbols
    tmp <- tmp %>%
      filter(rowSums(across(starts_with("p.value"), ~ . < 0.01), na.rm = TRUE) > 0) %>% 
      dplyr::mutate(across(contains("p.value"), ~ ifelse(is.na(.), "", gtools::stars.pval(.)))) |> 
      dplyr::mutate(across(contains("p.value"), ~ gsub("[.]", "", .)))
    
  } else if (FDR == T){
    
    # overall model FDR correct
    tmp1 <- lmm_data_overall %>%
      dplyr::filter(purrr::map_lgl(x, ~ any(stringr::str_detect(.x, exp)))) %>%
    
    tmp1_p <- tmp1 %>%
      dplyr::select(x,contains("p.value")) %>%
      pivot_longer(cols = contains("p.value"),names_to = "estimate", values_to = "p_value") %>%
      FDR_correct() %>%
      pivot_wider(id_cols = x, names_from = estimate, values_from = p_value)
    
    tmp1 <- tmp1 %>%
      select(!contains("p.value")) %>%
      left_join(tmp1_p, by = "x") %>%
      dplyr::select(x|contains(c("mpstatrsyes"))) %>% 
      dplyr::select(!contains(c("std.")))
    
    # interaction FDR correct
    tmp3 <- lmm_data_int %>%
      dplyr::filter(purrr::map_lgl(x, ~ any(stringr::str_detect(.x, exp)))) %>%
    
    tmp3_p <- tmp3 %>%
      dplyr::select(x,contains("p.value")) %>%
      pivot_longer(cols = contains("p.value"),names_to = "estimate", values_to = "p_value") %>%
      FDR_correct() %>%
      pivot_wider(id_cols = x, names_from = estimate, values_from = p_value)
    
    tmp3 <- tmp3 %>%
      select(!contains("p.value")) %>%
      left_join(tmp3_p, by = "x") %>%
      dplyr::select(!contains(c("std.", "age_at_consent", "bmi_at_consent", "interventionIdK",
                                "estimate_mpstatrsyes",
                                "p.value_mpstatrsyes"))) %>% 
      dplyr::select(!ends_with(c("visitIdM2", "visitIdM4", "visitIdM6")))
    
    # combine
    tmp <- tmp1 %>% 
      dplyr::full_join(tmp3, by = "x")
    
    #tmp$x = gsub(paste0(exp,"_"),"",tmp$x)
    
    # filter out only significant rows and replace pval with * symbols
    tmp <- tmp %>%
      filter(rowSums(across(starts_with("p.value"), ~ . < 0.05)) > 0, na.rm = TRUE) %>% 
      dplyr::mutate(across(contains("p.value"), ~ ifelse(is.na(.), "", gtools::stars.pval(.)))) |> 
      dplyr::mutate(across(contains("p.value"), ~ gsub("[.]", "", .)))
  }
  
  ## Split feature names into label/ome, fix names, order by ome ##----
  # split feature names into label/ome
  tmp$ome <- gsub("_.*", "", tmp$x)
  tmp$x <- mapply(function(pattern, text) gsub(pattern, "", text), tmp$ome, tmp$x) # remove ome names from x
  tmp$x <- gsub("^.", "",tmp$x) # remove trailing "_"
  
  # fix names clinical variables
  map <- read.csv("src/clinical_variables.csv")%>%
    dplyr::select(x, label) %>%
    dplyr::rename(feature_new = label)
  tmp <- merge(tmp, map, by = "x", all.x = TRUE)
  tmp$x <- ifelse(!is.na(tmp$feature_new), tmp$feature_new, tmp$x)
  tmp <- tmp %>% dplyr::select(-feature_new)
  
  # fix names immune cell populations
  load("src/populations_names_annotated.Rdata")
  populations = populations %>%
    dplyr::select(name, `population name`) %>%
    dplyr::rename(feature_new = `population name`, x = name)
  tmp <- merge(tmp, populations, by = "x", all.x = TRUE)
  tmp$x <- ifelse(!is.na(tmp$feature_new), tmp$feature_new, tmp$x)
  tmp <- tmp %>% dplyr::select(-feature_new)
  
  # fix names methylation scores
  indices_cerv <- readxl::read_xlsx("src/indices.xlsx", sheet = 1)
  indices_buccal <- readxl::read_xlsx("src/indices.xlsx", sheet = 2)
  indices_bl <- readxl::read_xlsx("src/indices.xlsx", sheet = 3)
  map <- rbind(indices_cerv,indices_buccal,indices_bl) %>%
    dplyr::distinct() %>%
    dplyr::rename(feature_new = label)
  tmp <- merge(tmp, map, by = "x", all.x = TRUE)
  tmp$x <- ifelse(!is.na(tmp$feature_new), tmp$feature_new, tmp$x)
  tmp <- tmp %>% dplyr::select(-feature_new)
  
  # order by ome
  tmp$ome <- gsub(": normalized", "", tmp$ome)
  tmp$ome <- factor(tmp$ome, levels = c("Blood haemogram","Body composition",
                                        "Flow cytometry: T cell staining","Flow cytometry: white blood cell staining",
                                        "Composite methylation scores: blood","Composite methylation scores: buccal","Composite methylation scores: cervical","Immune age: general",
                                        "Saliva nuclear magnetic resonance","Urine nuclear magnetic resonance",
                                        "Saliva microbiome: families", "Stool microbiome: families"))
  tmp <- tmp %>% 
    droplevels() %>%
    arrange(ome)
  
  ## Column splits and labels ##---- 
  col_split <- c("mpstatrs", "mpstatrs:time","mpstatrs:time","mpstatrs:time")
  col_split <- factor(col_split, levels = c("mpstatrs", "mpstatrs:time"))
  col_labs <- c("","M2", "M4", "M6")
  
  row_labs <- ifelse(tmp$`p.value_visitIdM6:mpstatrsyes` == "*" |
                     tmp$`p.value_visitIdM6:mpstatrsyes` == "**" |
                     tmp$`p.value_visitIdM6:mpstatrsyes` == "***" |
                     tmp$`p.value_visitIdM6:mpstatrsyes` == ".",
                     paste0("**", tmp$x, "**"),
                     tmp$x)
  
  # ome annotation
  color_palette <- carto_pal(length(levels(tmp$ome)), "Earth")
  color_omes <- setNames(color_palette, levels(tmp$ome))
  ha <- rowAnnotation(ome = tmp$ome,
                      col = list(ome = color_omes))
  
  ## Draw heatmap ##---- 
  
  # Indices
  ind_est <- grep("estimate", colnames(tmp))  
  ind_p <- grep("p.value", colnames(tmp))

  p <- Heatmap(as.matrix(tmp[,ind_est]),
               name = 'estimate\n(scaled)',
                
               # Row details
               row_labels = gt_render(row_labs),
               row_names_side = 'left',
                cluster_rows = F,
                cluster_row_slices = F,
               show_row_dend = F, 
               row_title = NULL,
               
                # Row annotation
               left_annotation = ha,
               
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
  
  
  return(p)
  
}
