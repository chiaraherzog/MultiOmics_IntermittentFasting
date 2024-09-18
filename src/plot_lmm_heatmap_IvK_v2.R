#' Plot a LMM Heatmap with Custom Annotations
#' 
#' @param lmm_data_int Data frame for the intervention assessment model.
#' @param cols Character vector of colors for the heatmap. Default colors are provided but can be customized.
#' @param relabel Data frame for relabeling variables in the plot. Must contain columns `x` (original variable) and `label` (new variable name).
#' @param filter_relabel Logical. If `TRUE`, only variables in the relabel data frame will be plotted. Default is `TRUE`.
#' @param relabel_assay Logical. If `TRUE`, assay names will be renamed based on the `assay2` column in the relabel data frame. Default is `FALSE`.
#' @param colour_assays List of colors for assays. Default colors are provided for clinical variables.
#' @param methyl_sub Character. Filter for methylation data from a specific tissue (e.g., "cervical", "buccal", "blood"). Default is `NULL`.
#' @param m6_sep Logical. If `TRUE`, M6 values will be plotted separately to avoid blank spots. Default is `TRUE`.
#' @param age_cor Logical. If `TRUE`, correlation with age will be plotted as an annotation. Default is `FALSE`.
#' @param cluster Character or logical. Determines row clustering method. Default is 'default', which orders by M6 time estimate. If `FALSE`, rows are not clustered.
#' @param pval_threshold Any row with an uncorrected p-value below threshold will be plotted. Default is 1.
#' @param top_n Within each assay, this number of significant top features will be kept. Default is 5.
#' @param pval_permutation_filter Logical. If `TRUE`, filters p-values using permutation tests. Default is `TRUE`.
#' @return A `ComplexHeatmap` object representing the LMM heatmap.
#' @import ComplexHeatmap
#' @import dplyr
#' @import tidyr
#' @import gtools
#' @import circlize
#' @export
plot_lmm_heatmap_IvK_v2 <- function(lmm_data_int,
                                cols = c("#1b69a1", "#48a0af", "#f39668", "#ec6669"),
                                relabel = NULL,
                                colour_assays = list(assay = c("Blood haemogram" = "#ec6669",
                                                               "Body weight and composition" = "#832c9b",
                                                               "Spirometry" = "#bd647d",
                                                               "Functional exercise capacity" = "#f39668",
                                                               "Blood test" = '#71b5a9',
                                                               "Skin histology and TEWL" = "#5f70a8",
                                                               "Vascular features" = '#0d49a1')),
                                relabel_assay = FALSE,
                                m6_sep = TRUE,
                                cluster = 'default',
                                pval_threshold = 1,
                                top_n = 5,
                                pval_permutation_filter = TRUE) {
  
  # Check for and install required packages
  if(!require("gtools")) install.packages("gtools")
  
  # Step 1: Filter relevant data frames and p-values, and apply permutation filtering if applicable #----
  tmp <- lmm_data_int |> 
    dplyr::select(!contains(c("std.", "age_at_consent", "bmi_at_consent", "compliancemedium",
                              "estimate_interventionIdK",
                              "p.value_interventionIdK"))) |> 
    dplyr::select(!ends_with(c("visitIdM2", "visitIdM4", "visitIdM6")))
  
  if(pval_permutation_filter) {
    source(here("src/getPermutationComparisonPval.R"))
    var_filter <- tmp$x[grepl(paste0(relabel$x, collapse = '|'), tmp$x)]
    
    load(here("out/outPermIntervention.Rdata"))
    tmp <- getPermutationComparisonPval(observed = tmp, random = outPermIntervent, variables = var_filter)
  }
  
  tmp <- tmp |> 
    dplyr::mutate(x = gsub("_clr_","_",x)) %>%
    tidyr::separate(x, into = c('assay', "x"), sep = "_", extra = 'merge')
  
  # Step2: Relabel variables if specified #----
  if(!is.null(relabel)) {
    tmp <- tmp |> 
      dplyr::left_join(relabel, by = c("x", "assay")) |> 
      dplyr::mutate(x = label) |> 
      dplyr::select(-label)
    
    if(relabel_assay == T) {
      tmp <- tmp |> dplyr::mutate(assay = assay2) |> dplyr::select(-assay2)
    }
  }
  
  # Step 3: Apply FDR correction, filter out only significant p-values based on p-value threshold, keeping only top n features within each assay, replace p-values by stars and bold significant p-values  #----
  source(here::here("src/FDRcorr.R"))
  tmpfdr <- FDRcorr(tmp, append = FALSE)
  
  keep = tmp %>%
    dplyr::filter_at(vars(contains("p.value")), any_vars(. < pval_threshold)) |> 
    dplyr::mutate(min_pvalue = do.call(pmin, across(contains("p.value")))) |> 
    dplyr::arrange(min_pvalue) |> 
    dplyr::group_by(assay) |> 
    dplyr::slice_head(n = top_n) |> 
    dplyr::ungroup() |>
    dplyr::pull(x2)
  
  p_value_columns <- grep("^p\\.value", names(tmp), value = TRUE)
  
  custom_stars_pval <- function(pvals) {
    sapply(pvals, function(pval) {
      if (is.na(pval)) {
        return("")  # Return an empty string for NA values
      } else {
        # Convert p-values to significance stars
        result <- cut(pval, breaks = c(0, 0.001, 0.01, 0.05, 1), 
                      labels = c("***", "**", "*", ""), right = FALSE)
        # Ensure the result is returned as a character
        return(as.character(result))
      }
    })
  }
  
  for (col in p_value_columns) {
    tmp[[col]] <- ifelse(is.na(tmpfdr[[col]]), "",
                         ifelse(tmpfdr[[col]] < 0.05,
                                paste0("<b>", custom_stars_pval(tmp[[col]]), "</b>"), 
                                custom_stars_pval(tmp[[col]])))
  }
  
  tmp = tmp %>% dplyr::filter(x2 %in% keep)

  # Set up indices and labels
  ind_est <- grep("estimate", colnames(tmp))  
  ind_p <- grep("p.value", colnames(tmp))
  
  col_labs <- c("M2", "M4", "M6")
  
  if(m6_sep) {
    row_split <- ifelse(is.na(tmp$`estimate_visitIdM2:interventionIdK`), 2, 1)
  } else {
    row_split <- rep(1, nrow(tmp))
  }
  
  row_labs <- ifelse(grepl("<b>", tmp$`p.value_visitIdM6:interventionIdK`, fixed = TRUE),
                     paste0("**", tmp$x, "**"),  # FDR significant (bold)
                     ifelse(grepl("*", tmp$`p.value_visitIdM6:interventionIdK`, fixed = TRUE),
                            paste0("*", tmp$x, "*"),  # Permutation significant (italic)
                            tmp$x))
  
  
  
  
  ha <- if(is.list(colour_assays)) {
    rowAnnotation(assay = tmp$assay, col = colour_assays)
  } else {
    NULL
  }
  
  hr <- NULL
  
  cell_fun <- function(j, i, x, y, width, height, fill) {
    full_text <- as.matrix(tmp[, ind_p])[i, j]
    if (grepl("<b>", full_text)) {
      clean_text <- gsub("<b>|</b>", "", full_text)
      grid.text(clean_text, x = x, y = y, 
                gp = gpar(fontsize = 7, fontface = "bold", fontfamily = "Arial Black"), just = "centre")
    } else {
      grid.text(full_text, x = x, y = y, 
                gp = gpar(fontsize = 7, fontface = "plain", col = 'grey20'), just = "centre")
    }
  }
  
  # Heatmap generation based on clustering option
  p <- Heatmap(
    as.matrix(tmp[, ind_est]),
    name = 'estimate\n(scaled)',
    row_labels = gt_render(row_labs),
    row_names_side = 'left',
    row_split = row_split,
    row_order = if(cluster == 'default') order(tmp$`p.value_visitIdM6:interventionIdK`, decreasing = TRUE) else NULL,
    cluster_rows = ifelse(cluster == FALSE, FALSE, TRUE),
    cluster_row_slices = FALSE,
    show_row_dend = FALSE,
    row_title = NULL,
    left_annotation = ha,
    right_annotation = hr,
    column_labels = col_labs,
    column_title_gp = grid::gpar(fontsize = 10),
    # column_split = col_split,
    cluster_columns = FALSE,
    show_column_dend = FALSE,
    cell_fun = cell_fun,
    na_col = 'white',
    col = circlize::colorRamp2(breaks = seq(-max(abs(na.omit(tmp[, ind_est]))),
                                            max(abs(na.omit(tmp[, ind_est]))),
                                            length.out = 30),
                               colors = colorRampPalette(c(cols[c(1, 2)], "grey95", cols[c(3, 4)]))(30)),
    row_names_gp = grid::gpar(fontsize = 9),
    column_names_gp = grid::gpar(fontsize = 9),
    border_gp = gpar(lwd = 0.5),
    border = TRUE
  )
  
  # # Add age correlation line if applicable
  # if(age_cor && cluster != FALSE) {
  #   p <- p + decorate_annotation("Correlation with age\n(baseline)", {
  #     grid.lines(c(0.5), gp = gpar(col = "grey80", lty = 'dotted'))
  #   })
  # }
  
  return(p)
}
