<<<<<<< HEAD
#' @param lmm_data_time Minimal model that provides time coefficient (no interaction) (out_lmm$`Minimal model`)
#' @param lmm_data_compliance Basic model that provides time*compliance interaction (out_lmm$`Basic model with BMI`)
#' @param cols Colour palette (leave blank for default colours [recommended], but can be altered if desired)
#' @param relabel (optional). Data frame supplying labels for the plot. x = variable, label = relabelled variable
#' @param filter_relabel (optional) boolean: if supplying new labels, should only variables provided in that dataframe be plotted? default to TRUE
#' @param relabel_assay (optional). Set to TRUE if you would like to rename the assay. Assay names should be provided in the relabel parameter, under column name 'assay2'
#' @param colour_assays (optional) if not NULL, list of default colours will be used for assays. This only works for clinical variables so far.
#' @param methyl_sub default to NULL; filter for methylation data of one tissue only (cervical, buccal, or blood)
#' @return ComplexHeatmap object


plot_lmm_heatmap_v5 <- function(lmm_data_time,
                                lmm_data_compliance,
                                cols = c("#1b69a1",
                                         "#48a0af",
                                         "#f39668",
                                         "#ec6669"),
                                relabel = NULL,
                                filter_relabel = T,
                                colour_assays = NULL,
                                relabel_assay = NULL,
                                methyl_sub = NULL){
  
  if(!require("gtools")){
    install.packages("gtools")
  }
  
  tmp1 <- lmm_data_time |> 
    dplyr::select(!contains(c("std.", "age_at_consent", "bmi_at_consent"))) |>
    # set pval
    dplyr::mutate(across(contains("p.value"), ~ ifelse(is.na(.), " ", gtools::stars.pval(.)))) |> 
    dplyr::mutate(across(contains("p.value"), ~ gsub("[.]", " ", .)))
  
  tmp2 <- lmm_data_compliance |> 
    dplyr::select(!contains(c("std", "compliancemedium", "bmi_at_consent", 
                              "age_at_consent", "estimate_compliancehigh", "p.value_compliancehigh"))) |> 
    dplyr::select(!ends_with(c("visitIdM2", "visitIdM4", "visitIdM6"))) |> 
    
    # set pval
    dplyr::mutate(across(contains("p.value"), ~ ifelse(is.na(.), " ", gtools::stars.pval(.)))) |> 
    dplyr::mutate(across(contains("p.value"), ~ gsub("[.]", " ", .)))
  
  tmp <- tmp1 |> 
    dplyr::left_join(tmp2)
  
  # condition: if assay name in x, need to split it for relabelling to work
  if(any(grepl("Composite|Functional|Blood haemogram|Flow cytometry|Body composition", tmp$x))){
    
    tmp <- tmp |> 
      tidyr::separate(x, into = c('assay', "x"),
                      sep = "_", 
                      extra = 'merge')
    
    # methyl subset
    if(!is.null(methyl_sub)){
      tmp <- tmp |> 
        dplyr::filter(grepl(methyl_sub, assay) & grepl("methylation", assay))
    }
    
  }
  
  if(!is.null(relabel)){
    
    if(filter_relabel == T){
      tmp <- tmp |> 
        dplyr::filter(x %in% relabel$x)
      
      relabel <- relabel[relabel$x %in% tmp$x,]
      
      tmp <- tmp[match(relabel$x, tmp$x),]
    }
    
    tmp <- tmp |> 
      dplyr::left_join(relabel) |> 
      dplyr::mutate(x = label) |> 
      dplyr::select(-label)
    
    if(relabel_assay == T){
      tmp <- tmp |> 
        dplyr::mutate(assay = assay2)
    }
    
  }
  
  # reorder rows by assay order provided in colour_assays
  if("assay" %in% colnames(tmp) & is.list(colour_assays)){
    
    assay_order <- names(colour_assays$assay)
    
    tmp <- tmp |> 
      dplyr::mutate(
        assay = factor(assay,
                       levels = assay_order)
      ) |> 
      dplyr::arrange(assay, x)
  }
  
  # Indices
  ind_est <- grep("estimate", colnames(tmp))  
  ind_p <- grep("p.value", colnames(tmp))
  
  # Column splits and labels
  col_split <- factor(rep(c("time", "high\ncompliance"), each = 3), levels = c("time", "high\ncompliance"))
  col_labs <- rep(c("M2", "M4", "M6"), 2)
  
  row_labs <- tmp$x
  
  if(is.list(colour_assays)){
    ha <- rowAnnotation(
      type = tmp$assay,
      col = list(type = colour_assays$assay)
    )
  } else {
    ha <- NULL
  }
  
  hr = NULL
  
  p <- Heatmap(
    as.matrix(tmp[, ind_est]),
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
    right_annotation = hr,
    
    # Column details
    column_labels = col_labs,
    column_title_gp = grid::gpar(fontsize = 10),
    column_split = col_split,
    cluster_columns = F,
    show_column_dend = F,
    
    # Annotation
    cell_fun = function(j, i, x, y, width, height, fill) {
      grid.text(as.matrix(tmp[, ind_p])[i, j],
                x,
                y,
                gp = gpar(fontsize = 7),
                just = "centre")
    },
    
    # Colours
    na_col = 'white',
    col = circlize::colorRamp2(
      breaks = seq(-max(abs(na.omit(
        tmp[, ind_est]
      ))), max(abs(na.omit(
        tmp[, ind_est]
      ))), length.out = 30),
      colors = colorRampPalette(c(cols[c(1, 2)], "grey95", cols[c(3, 4)]))(30)
      #                            # colors = rev(color("batlow")(5)),
      #                            # colors = rev(cols[c(2,5,4, 6, 7)])
    ),
    
    # Titles
    row_names_gp = grid::gpar(fontsize = 9),
    column_names_gp = grid::gpar(fontsize = 9),
    
    # Borders
    border_gp = gpar(lwd = 0.5),
    border = T
  )
  
  full_df <- cbind(tmp$x,
                   tmp[,ind_est], tmp[,ind_p])
  
  out <- list(val = tmp[,ind_est],
              pval = tmp[,ind_p],
              plot = p,
              full_df = full_df)
  
  
  return(out)
  
}
=======
#' @param lmm_data_time Minimal model that provides time coefficient (no interaction) (out_lmm$`Minimal model`)
#' @param lmm_data_compliance Basic model that provides time*compliance interaction (out_lmm$`Basic model with BMI`)
#' @param cols Colour palette (leave blank for default colours [recommended], but can be altered if desired)
#' @param relabel (optional). Data frame supplying labels for the plot. x = variable, label = relabelled variable
#' @param filter_relabel (optional) boolean: if supplying new labels, should only variables provided in that dataframe be plotted? default to TRUE
#' @param relabel_assay (optional). Set to TRUE if you would like to rename the assay. Assay names should be provided in the relabel parameter, under column name 'assay2'
#' @param colour_assays (optional) if not NULL, list of default colours will be used for assays. This only works for clinical variables so far.
#' @param methyl_sub default to NULL; filter for methylation data of one tissue only (cervical, buccal, or blood)
#' @return ComplexHeatmap object


plot_lmm_heatmap_v5 <- function(lmm_data_time,
                                lmm_data_compliance,
                                cols = c("#1b69a1",
                                         "#48a0af",
                                         "#f39668",
                                         "#ec6669"),
                                relabel = NULL,
                                filter_relabel = T,
                                colour_assays = NULL,
                                relabel_assay = NULL,
                                methyl_sub = NULL){
  
  if(!require("gtools")){
    install.packages("gtools")
  }
  
  tmp1 <- lmm_data_time |> 
    dplyr::select(!contains(c("std.", "age_at_consent", "bmi_at_consent"))) |>
    # set pval
    dplyr::mutate(across(contains("p.value"), ~ ifelse(is.na(.), " ", gtools::stars.pval(.)))) |> 
    dplyr::mutate(across(contains("p.value"), ~ gsub("[.]", " ", .)))
  
  tmp2 <- lmm_data_compliance |> 
    dplyr::select(!contains(c("std", "compliancemedium", "bmi_at_consent", 
                              "age_at_consent", "estimate_compliancehigh", "p.value_compliancehigh"))) |> 
    dplyr::select(!ends_with(c("visitIdM2", "visitIdM4", "visitIdM6"))) |> 
    
    # set pval
    dplyr::mutate(across(contains("p.value"), ~ ifelse(is.na(.), " ", gtools::stars.pval(.)))) |> 
    dplyr::mutate(across(contains("p.value"), ~ gsub("[.]", " ", .)))
  
  tmp <- tmp1 |> 
    dplyr::left_join(tmp2)
  
  # condition: if assay name in x, need to split it for relabelling to work
  if(any(grepl("Composite|Functional|Blood haemogram|Flow cytometry|Body composition", tmp$x))){
    
    tmp <- tmp |> 
      tidyr::separate(x, into = c('assay', "x"),
                      sep = "_", 
                      extra = 'merge')
    
    # methyl subset
    if(!is.null(methyl_sub)){
      tmp <- tmp |> 
        dplyr::filter(grepl(methyl_sub, assay) & grepl("methylation", assay))
    }
    
  }
  
  if(!is.null(relabel)){
    
    if(filter_relabel == T){
      tmp <- tmp |> 
        dplyr::filter(x %in% relabel$x)
      
      relabel <- relabel[relabel$x %in% tmp$x,]
      
      tmp <- tmp[match(relabel$x, tmp$x),]
    }
    
    tmp <- tmp |> 
      dplyr::left_join(relabel) |> 
      dplyr::mutate(x = label) |> 
      dplyr::select(-label)
    
    if(relabel_assay == T){
      tmp <- tmp |> 
        dplyr::mutate(assay = assay2)
    }
    
  }
  
  # reorder rows by assay order provided in colour_assays
  if("assay" %in% colnames(tmp) & is.list(colour_assays)){
    
    assay_order <- names(colour_assays$assay)
    
    tmp <- tmp |> 
      dplyr::mutate(
        assay = factor(assay,
                       levels = assay_order)
      ) |> 
      dplyr::arrange(assay, x)
  }
  
  # Indices
  ind_est <- grep("estimate", colnames(tmp))  
  ind_p <- grep("p.value", colnames(tmp))
  
  # Column splits and labels
  col_split <- factor(rep(c("time", "high\ncompliance"), each = 3), levels = c("time", "high\ncompliance"))
  col_labs <- rep(c("M2", "M4", "M6"), 2)
  
  row_labs <- tmp$x
  
  if(is.list(colour_assays)){
    ha <- rowAnnotation(
      type = tmp$assay,
      col = list(type = colour_assays$assay)
    )
  } else {
    ha <- NULL
  }
  
  hr = NULL
  
  p <- Heatmap(
    as.matrix(tmp[, ind_est]),
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
    right_annotation = hr,
    
    # Column details
    column_labels = col_labs,
    column_title_gp = grid::gpar(fontsize = 10),
    column_split = col_split,
    cluster_columns = F,
    show_column_dend = F,
    
    # Annotation
    cell_fun = function(j, i, x, y, width, height, fill) {
      grid.text(as.matrix(tmp[, ind_p])[i, j],
                x,
                y,
                gp = gpar(fontsize = 7),
                just = "centre")
    },
    
    # Colours
    na_col = 'white',
    col = circlize::colorRamp2(
      breaks = seq(-max(abs(na.omit(
        tmp[, ind_est]
      ))), max(abs(na.omit(
        tmp[, ind_est]
      ))), length.out = 30),
      colors = colorRampPalette(c(cols[c(1, 2)], "grey95", cols[c(3, 4)]))(30)
      #                            # colors = rev(color("batlow")(5)),
      #                            # colors = rev(cols[c(2,5,4, 6, 7)])
    ),
    
    # Titles
    row_names_gp = grid::gpar(fontsize = 9),
    column_names_gp = grid::gpar(fontsize = 9),
    
    # Borders
    border_gp = gpar(lwd = 0.5),
    border = T
  )
  
  full_df <- cbind(tmp$x,
                   tmp[,ind_est], tmp[,ind_p])
  
  out <- list(val = tmp[,ind_est],
              pval = tmp[,ind_p],
              plot = p,
              full_df = full_df)
  
  
  return(out)
  
}
>>>>>>> ccb3274a4298533190659db2ca2fb40870057cca
