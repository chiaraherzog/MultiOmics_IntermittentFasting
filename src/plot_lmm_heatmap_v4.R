#' Plot a LMM Heatmap with Custom Annotations
#' 
#' @param lmm_data_time Data frame for the minimal model that provides time coefficients (no interaction).
#' @param lmm_data_compliance Data frame for the basic model that provides time*compliance interaction.
#' @param cols Character vector of colors for the heatmap. Default colors are provided but can be customized.
#' @param relabel Data frame for relabeling variables in the plot. Must contain columns `x` (original variable) and `label` (new variable name).
#' @param filter_relabel Logical. If `TRUE`, only variables in the relabel data frame will be plotted. Default is `TRUE`.
#' @param relabel_assay Logical. If `TRUE`, assay names will be renamed based on the `assay2` column in the relabel data frame. Default is `FALSE`.
#' @param colour_assays List of colors for assays. Default colors are provided for clinical variables.
#' @param methyl_sub Character. Filter for methylation data from a specific tissue (e.g., "cervical", "buccal", "blood"). Default is `NULL`.
#' @param m6_sep Logical. If `TRUE`, M6 values will be plotted separately to avoid blank spots. Default is `TRUE`.
#' @param age_cor Logical. If `TRUE`, correlation with age will be plotted as an annotation. Default is `FALSE`.
#' @param cluster Character or logical. Determines row clustering method. Default is 'default', which orders by M6 time estimate. If `FALSE`, rows are not clustered.
#' @param pval_permutation_filter Logical. If `TRUE`, filters p-values using permutation tests. Default is `TRUE`.
#' @return A `ComplexHeatmap` object representing the LMM heatmap.
#' @import ComplexHeatmap
#' @import dplyr
#' @import tidyr
#' @import gtools
#' @import circlize
#' @export
plot_lmm_heatmap_v4 <- function(lmm_data_time,
                                lmm_data_compliance,
                                cols = c("#1b69a1", "#48a0af", "#f39668", "#ec6669"),
                                relabel = NULL,
                                filter_relabel = TRUE,
                                colour_assays = list(assay = c("Blood haemogram" = "#ec6669",
                                                               "Body weight and composition" = "#832c9b",
                                                               "Spirometry" = "#bd647d",
                                                               "Functional exercise capacity" = "#f39668",
                                                               "Blood test" = '#71b5a9',
                                                               "Skin histology and TEWL" = "#5f70a8",
                                                               "Vascular features" = '#0d49a1')),
                                relabel_assay = FALSE,
                                m6_sep = TRUE,
                                methyl_sub = NULL,
                                age_cor = FALSE,
                                cluster = 'default',
                                pval_permutation_filter = TRUE) {
  
  # Check for and install required packages
  if(!require("gtools")) install.packages("gtools")
  
  # Step 1: Filter relevant data frames and p-values, and apply permutation filtering if applicable
  tmp1 <- lmm_data_time |> 
    dplyr::select(!contains(c("std.", "age_at_consent", "bmi_at_consent")))
  
  if(pval_permutation_filter) {
    source(here("src/getPermutationComparisonPval.R"))
    var_filter <- tmp1$x[grepl(paste0(relabel$x, collapse = '|'), tmp1$x)]
    
    load(here("out/outPermMinimal.Rdata"))
    tmp1 <- getPermutationComparisonPval(observed = tmp1, random = outPermMinimal, variables = var_filter)
  }
  
  tmp2 <- lmm_data_compliance |> 
    dplyr::select(!contains(c("std", "compliancemedium", "bmi_at_consent", 
                              "age_at_consent", "estimate_compliancehigh", "p.value_compliancehigh"))) |> 
    dplyr::select(!ends_with(c("visitIdM2", "visitIdM4", "visitIdM6")))
  
  if(pval_permutation_filter) {
    load("out/outPerm.Rdata")
    tmp2 <- getPermutationComparisonPval(observed = tmp2, random = outPerm, variables = var_filter)
  }
  
  # Merge time and interaction data
  tmp <- tmp1 |> 
    dplyr::left_join(tmp2) |> 
    tidyr::separate(x, into = c('assay', "x"), sep = "_", extra = 'merge')
  
  # Filter for methylation data if specified
  if(!is.null(methyl_sub)) {
    tmp <- tmp |> dplyr::filter(grepl(methyl_sub, assay) & grepl("methylation", assay))
  }
  
  # Relabel variables if specified
  if(!is.null(relabel)) {
    if(filter_relabel) {
      if(any(grepl("stimul", relabel$assay))){
        tmp <- tmp |>
          dplyr::inner_join(relabel)
      } else {
      intersect <- intersect(tmp$x, relabel$x)
      relabel <- relabel[relabel$x %in% intersect,]
      tmp <- tmp[match(relabel$x, tmp$x),]
      }
    }
  }
    
    if(age_cor == T) {
      load(here("out/corrAgeBaseline.R"))
      tmp$fullname <- paste0(tmp$assay,"_", tmp$x)
      corr <- corr |> 
        tidyr::separate(rowname, "_", into = c("assay", "name"), extra = 'merge',remove = F) |> 
        dplyr::slice(match(tmp$fullname, rowname))
      
      tmp <- tmp |>
        dplyr::left_join(dplyr::select(corr, rowname, name, cor, p), by = c('x' = 'name',
                                                                                         'fullname' = 'rowname')) |> 
        dplyr::select(-c(fullname))
      
    }
    
    tmp <- tmp |> 
      dplyr::left_join(relabel) |> 
      dplyr::mutate(x = label) |> 
      dplyr::select(-label)
    
    if(relabel_assay == T) {
      tmp <- tmp |> dplyr::mutate(assay = assay2) |> dplyr::select(-assay2)
    }
  
  # Apply FDR correction and bold significant p-values
  source(here::here("src/FDRcorr.R"))
  tmpfdr <- FDRcorr(tmp, append = FALSE)
  
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
  
  
  # Set up indices and labels
  ind_est <- grep("estimate", colnames(tmp))  
  ind_p <- grep("p.value", colnames(tmp))
  
  col_split <- factor(rep(c("time", "high\ncompliance"), each = 3), levels = c("time", "high\ncompliance"))
  col_labs <- rep(c("M2", "M4", "M6"), 2)
  
  if(m6_sep) {
    row_split <- ifelse(is.na(tmp$`estimate_visitIdM2`), 2, 1)
  } else {
    row_split <- rep(1, nrow(tmp))
  }
  
  row_labs <- ifelse(grepl("<b>", tmp$`p.value_visitIdM6:compliancehigh`, fixed = TRUE) | 
                       grepl("<b>", tmp$p.value_visitIdM6, fixed = TRUE),
                     paste0("**", tmp$x, "**"),  # FDR significant (bold)
                     ifelse(grepl("*", tmp$`p.value_visitIdM6:compliancehigh`, fixed = TRUE) | 
                              grepl("*", tmp$p.value_visitIdM6, fixed = TRUE),
                            paste0("*", tmp$x, "*"),  # Permutation significant (italic)
                            tmp$x))
  
  
  
  
  ha <- if(is.list(colour_assays)) {
    rowAnnotation(assay = tmp$assay, col = colour_assays)
  } else {
    NULL
  }
  
  hr <- if(age_cor) {
    rowAnnotation(
      'Correlation with age\n(baseline)' = anno_points(tmp$cor,
                                                       pch = ifelse(tmp$p < 0.05, 19, 1),
                                                       ylim = c(-1, 1),
                                                       gp = gpar(col = ifelse(tmp$cor < 0, cols[1], cols[4]))),
      'Opposite change\nwith\nintervention' = anno_text(ifelse(sign(tmp$cor) != sign(tmp$estimate_visitIdM6) & tmp$p.value_visitIdM6 != " ", "←", ""))
    )
  } else {
    NULL
  }
  
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
    row_order = if(cluster == 'default') order(tmp$estimate_visitIdM6, decreasing = TRUE) else NULL,
    cluster_rows = ifelse(cluster == FALSE, FALSE, TRUE),
    cluster_row_slices = FALSE,
    show_row_dend = FALSE,
    row_title = NULL,
    left_annotation = ha,
    right_annotation = hr,
    column_labels = col_labs,
    column_title_gp = grid::gpar(fontsize = 10),
    column_split = col_split,
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
