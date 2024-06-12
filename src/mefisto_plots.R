
# custom plots for mefisto trained models
# version for IF
# time/compliance box-plots depends on paired_longitudinal_compliance.R
# author: Charlotte Vavourakis

cols <- c("#1b69a1", "#48a0af", "#71b5a9", "#ec6669", "#f39668", "#bd647d", "#832c9b", "#5f70a8")

col_visitId <- c("M0" = cols[2], "M2" = cols[1], "M4" = cols[4], "M6" = cols[5])

# PCA-like plot, first two factors, color by time #----

mefisto_biplot_F12 <- function(df1, df2, time_type, time_unit = NULL) {
  
  library(tidyverse)
  library(rcartocolor)
  
  # rename groups
  
  df1 = df1 %>%
    dplyr::mutate(group = case_when(
      group == "I" ~ "no MCT",
      group == "K" ~ "+ MCT"))
  
  df2 = df2 %>%
    dplyr::mutate(group = case_when(
      group == "I" ~ "no MCT",
      group == "K" ~ "+ MCT"))
  
  if (time_type == "factor") {
    tmp = df2 %>%
      dplyr::filter(factor %in% c("Factor1","Factor2")) %>%
      dplyr::select(sample,factor,value.factor,color_by, group) %>%
      dplyr::rename(time=color_by) %>%
      tidyr::pivot_wider(names_from = factor, values_from = value.factor)
    
    xl = df1 %>%
      dplyr::filter(factor == "Factor1") %>%
      pull(value) %>%
      mean() %>%
      round(digits = 1)
    yl = df1 %>%
      dplyr::filter(factor == "Factor2") %>%
      pull(value) %>%
      mean() %>%
      round(digits = 1)
    
    if (length(unique(tmp$group)) == 1) {
      plot = tmp %>%
        ggplot(aes(x = Factor1, 
                   y = Factor2, 
                   fill = time, 
                   color = time)) +
        geom_point(shape = 21, alpha = 0.5) + 
        theme_minimal() +
        stat_ellipse(type="norm", show.legend = FALSE) +
        scale_fill_manual(values = col_visitId) +
        scale_color_manual(values = col_visitId, guide = "none") +
        xlab(paste0("Factor 1\n(",xl,"% mean variance explained)")) +
        ylab(paste0("Factor 2\n(",yl,"% mean variance explained)")) +
        labs(fill="")
    } else {
      plot = tmp %>%
        ggplot(aes(x = Factor1, 
                   y = Factor2, 
                   fill = time, 
                   color = time)) +
        geom_point(shape = 21, alpha = 0.5) + 
        theme_bw() +
        stat_ellipse(type="norm", show.legend = FALSE) +
        scale_fill_manual(values = col_visitId) +
        scale_color_manual(values = col_visitId, guide = "none") +
        xlab(paste0("Factor 1\n(",xl,"% mean variance explained)")) +
        ylab(paste0("Factor 2\n(",yl,"% mean variance explained)")) +
        labs(fill="") +
        facet_wrap(~group) +
        theme(strip.background = element_blank(),  # Remove strip background
              strip.text = element_text(hjust = 0))
    }

  } else if (time_type == "continuous"){
    
    tmp = df2 %>%
      dplyr::filter(factor %in% c("Factor1","Factor2")) %>%
      dplyr::select(sample,factor,value.factor,value.covariate) %>%
      dplyr::rename(time=value.covariate) %>%
      tidyr::pivot_wider(names_from = factor, values_from = value.factor)
    
    xl = df1 %>%
      dplyr::filter(factor == "Factor1") %>%
      pull(value) %>%
      mean() %>%
      round(digits = 1)
    yl = df1 %>%
      dplyr::filter(factor == "Factor2") %>%
      pull(value) %>%
      mean() %>%
      round(digits = 1)
    
    
    plot = tmp %>%
      ggplot(aes(x = Factor1, 
                 y = Factor2, 
                 fill = time, 
                 color = time)) +
      geom_point(color = "black", shape = 21) + 
      theme_minimal() +
      scale_fill_carto_c(palette = "Earth") +
      xlab(paste0("Factor 1\n",xl,"% mean variance explained")) +
      ylab(paste0("Factor 2\n",yl,"% mean variance explained")) +
      labs(fill=paste0("time (",time_unit,")"))
  }
  
  
  
  return(plot)
}


# factor variance and significance association with selected covariates #----

mefisto_factor_covcor <- function(df1, df2, feat, labs){
  
  library(tidyverse)
  library(ComplexHeatmap)
  
  fs = df2 %>%
    tidyr::pivot_wider(names_from = factor, values_from = value.factor) %>%
    dplyr::select(c(sample,starts_with("Factor")))
  
  nfac = ncol(fs) -1
  fs = fs %>%
    dplyr::select(c(sample,paste0("Factor",seq(1:nfac))))
  
  tmp = metadat %>% 
    dplyr::select(sample,any_of(features)) %>% 
    dplyr::mutate_if(is.character, as.factor) %>%
    dplyr::mutate(mpstatrs = gsub("[*]", "", mpstatrs)) %>%
    dplyr::select(where(~n_distinct(.) >= 2)) %>%
    dplyr::full_join(fs, by = "sample")
  
  fs = tmp %>%
    dplyr::select(c(sample,starts_with("Factor")))
  
  tmp = tmp %>%
    dplyr::select(-starts_with("Factor"))
  
  identical(as.character(tmp$sample),as.character(fs$sample))
  
  fs = fs %>% dplyr::select(-sample)
  tmp = tmp %>% dplyr::select(-sample)
  
  mat = matrix(ncol = ncol(fs),
               nrow = ncol(tmp))
  colnames(mat) = colnames(fs)
  rownames(mat) = colnames(tmp)
  
  for (i in colnames(mat)){
    for (j in rownames(mat)){
      if(is.numeric(tmp[,j])){
        mat[j,i] <- cor.test(fs[,i], tmp[,j])$p.value
        # pmat[j,i] <- cor.test(pcs[,i], tmp[,j])
      } else {
        mat[j,i] <- kruskal.test(fs[,i], tmp[,j])$p.value
      }
    }
  }
  
  # for the other clinical data need to deal with NAs
  map = mat
  map = apply(map, 2, function(t) ifelse(t < 0.05, t, NA))
  map = map[,1:ncol(fs)]
  
  plot = Heatmap(-log10(map),
                 row_labels = labs,
                 cluster_columns = F,
                 row_names_side = 'left',
                 cluster_rows = F,
                 column_title = NULL,
                 show_row_dend = F, show_column_dend = F,
                 name = '-log10(p value)',
                 na_col = 'white',
                 col = circlize::colorRamp2(breaks = seq(40, 1.3, length.out = 5),
                                            colors = rev(viridis::viridis(5)),
                 ),
                 row_names_gp = grid::gpar(fontsize = 9),
                 column_names_gp = grid::gpar(fontsize = 9),
                 border_gp = gpar(lwd = 0.5),
                 border = T)
  
  return(plot)
  
}


# time/compliance box-plots values selected factor #----
mefisto_factor_box <- function(df2, metadat, factor_select, high_only = F) {
  
  library(tidyverse)
  library(patchwork)
  
  colours = cols[c(4,5,1)]
  
  if (high_only == F) {
    
    tmp = df2 %>%
      dplyr::filter(factor == factor_select) %>%
      dplyr::left_join(metadat, by = c("sample","group")) %>%
      dplyr::select(sample, subjectId, visitId, compliance, value.factor,group)
    
    tmp$compliance = factor(tmp$compliance, levels = c("low", "medium", "high"))
    
    complete = tmp |> 
      dplyr::group_by(subjectId) |> 
      dplyr::count() |> 
      ungroup() |> 
      dplyr::filter(n == 4)
    
    if (length(unique(tmp$group)) == 1) {
      
      plot = tmp %>%
        dplyr::filter(subjectId %in% c(complete$subjectId)) %>%
        dplyr::mutate(complabel = paste0(compliance, "<br>compliance"),
                      complabel = factor(complabel, levels = c("low<br>compliance", 'medium<br>compliance', 'high<br>compliance'))) %>% 
        ggplot(aes(x = visitId,
                   y = value.factor)) +
        geom_boxplot(alpha = 0.3,
                     aes(fill = compliance)) +
        geom_line(aes(group = subjectId,
                      colour = compliance),
                  alpha = 0.3) +
        ggpubr::stat_compare_means(ref.group = "M0",
                                   label = 'p.signif',
                                   size = 2.7,
                                   paired = T,
                                   label.y.npc = 0.95) +
        theme_bw() +
        theme(legend.position = 'none',
              aspect.ratio = 2,
              axis.title.y = ggtext::element_markdown(),
              strip.text.x = ggtext::element_markdown()) +
        labs(x = "", y = paste0(factor_select," value")) +
        scale_colour_manual(values = colours,
                            aesthetics = c("fill", 'colour')) +
        facet_wrap(~complabel) 
    } else {
      
      tmp1 = tmp %>% dplyr::filter(group == "I") 
      tmp2 = tmp %>% dplyr::filter(group == "K") 
      
      plot1 = tmp1 %>%
        dplyr::filter(subjectId %in% c(complete$subjectId)) %>%
        dplyr::mutate(complabel = paste0(compliance, "<br>compliance"),
                      complabel = factor(complabel, levels = c("low<br>compliance", 'medium<br>compliance', 'high<br>compliance'))) %>% 
        ggplot(aes(x = visitId,
                   y = value.factor)) +
        geom_boxplot(alpha = 0.3,
                     aes(fill = compliance)) +
        geom_line(aes(group = subjectId,
                      colour = compliance),
                  alpha = 0.3) +
        ggpubr::stat_compare_means(ref.group = "M0",
                                   label = 'p.signif',
                                   size = 2.7,
                                   paired = T,
                                   label.y.npc = 0.93) +
        theme_bw() +
        theme(legend.position = 'none',
              #aspect.ratio = 3,
              axis.title.y = ggtext::element_markdown(),
              strip.text.x = ggtext::element_markdown()) +
        labs(x = "", y = paste0(factor_select," value")) +
        scale_colour_manual(values = colours,
                            aesthetics = c("fill", 'colour')) +
        facet_wrap(~complabel) +
        ggtitle("no MCT")
      
      plot2 = tmp2 %>%
        dplyr::filter(subjectId %in% c(complete$subjectId)) %>%
        dplyr::mutate(complabel = paste0(compliance, "<br>compliance"),
                      complabel = factor(complabel, levels = c("low<br>compliance", 'medium<br>compliance', 'high<br>compliance'))) %>% 
        ggplot(aes(x = visitId,
                   y = value.factor)) +
        geom_boxplot(alpha = 0.3,
                     aes(fill = compliance)) +
        geom_line(aes(group = subjectId,
                      colour = compliance),
                  alpha = 0.3) +
        ggpubr::stat_compare_means(ref.group = "M0",
                                   label = 'p.signif',
                                   size = 2.7,
                                   paired = T,
                                   label.y.npc = 0.93) +
        theme_bw() +
        theme(legend.position = 'none',
              #aspect.ratio = 3,
              axis.title.y = ggtext::element_markdown(),
              strip.text.x = ggtext::element_markdown()) +
        labs(x = "", y = paste0(factor_select," value")) +
        scale_colour_manual(values = colours,
                            aesthetics = c("fill", 'colour')) +
        facet_wrap(~complabel) +
        ggtitle("+ MCT")
      
      plot = plot2 + plot1
    }
  } else {
    
    tmp = df2 %>%
      dplyr::filter(factor == factor_select) %>%
      dplyr::left_join(metadat, by = c("sample","group")) %>%
      dplyr::select(sample, subjectId, visitId, value.factor,group)
    
    complete = tmp |> 
      dplyr::group_by(subjectId) |> 
      dplyr::count() |> 
      ungroup() |> 
      dplyr::filter(n == 4)
  
      plot = tmp %>%
        dplyr::filter(subjectId %in% c(complete$subjectId)) %>%
        ggplot(aes(x = visitId,
                   y = value.factor)) +
        geom_boxplot(alpha = 0.3,
                     fill = colours[3]) +
        geom_line(aes(group = subjectId),
                  colour = colours[3],
                  alpha = 0.3) +
        ggpubr::stat_compare_means(ref.group = "M0",
                                   label = 'p.signif',
                                   size = 2.7,
                                   paired = T,
                                   label.y.npc = 0.95) +
        theme_bw() +
        theme(legend.position = 'none',
              aspect.ratio = 2,
              axis.title.y = ggtext::element_markdown(),
              strip.text.x = ggtext::element_markdown()) +
        labs(x = "", y = paste0(factor_select," value")) 
  }
  
  return(plot)
  
}
# coefficient plot top x weights selected factor #----

mefisto_plot_weights <- function(w_f, factor_select, n_feat, map_tcs_path, map_wb_path){
  
  library(tidyverse)
  library(openxlsx)
  library(rcartocolor)
  
  tmp = w_f %>%
    dplyr::arrange(desc(abs(value))) %>% 
    dplyr::slice(1:n_feat)
  
  # fix feature names
  tmp$label <- gsub(".*_", "", tmp$feature)
  tmp$feature <- gsub("_[^_]*$", "", tmp$feature)

  ## FCT data
  map_tcs = readRDS(map_tcs_path) %>% dplyr::rename(feature_new = feature,feature = featureId_mefisto)
  map_wb = readRDS(map_wb_path) %>% dplyr::rename(feature_new = feature,feature = featureId_mefisto)
  
  tmp <- merge(tmp, map_tcs, by = "feature", all.x = TRUE)
  tmp$feature <- ifelse(!is.na(tmp$feature_new), tmp$feature_new, tmp$feature)
  tmp <- tmp %>% dplyr::select(-feature_new)
  tmp <- merge(tmp, map_wb, by = "feature", all.x = TRUE)
  tmp$feature <- ifelse(!is.na(tmp$feature_new), tmp$feature_new, tmp$feature)
  tmp <- tmp %>% dplyr::select(-feature_new)
  
  load("src/populations_names_annotated.Rdata")
  populations = populations %>%
    dplyr::select(name, `population name`) %>%
    dplyr::rename(feature_new = `population name`, feature = name)
  tmp <- merge(tmp, populations, by = "feature", all.x = TRUE)
  tmp$feature <- ifelse(!is.na(tmp$feature_new), tmp$feature_new, tmp$feature)
  tmp <- tmp %>% dplyr::select(-feature_new)
  
  ## methylation scores
  indices_cerv <- readxl::read_xlsx("src/indices.xlsx", sheet = 1)
  indices_buccal <- readxl::read_xlsx("src/indices.xlsx", sheet = 2)
  indices_bl <- readxl::read_xlsx("src/indices.xlsx", sheet = 3)
  map <- rbind(indices_cerv,indices_buccal,indices_bl) %>%
    dplyr::distinct() %>%
    dplyr::rename(feature_new = label, feature = x)
  
  tmp <- merge(tmp, map, by = "feature", all.x = TRUE)
  tmp$feature <- ifelse(!is.na(tmp$feature_new), tmp$feature_new, tmp$feature)
  tmp <- tmp %>% dplyr::select(-feature_new)
  
  ## paste label back, avoiding ggplot complaining about duplicate labels
  tmp$feature <- paste0(tmp$feature,"_",tmp$label)
  
  # fix order
  tmp <- tmp %>% arrange(view)
  tmp$feature <- factor(tmp$feature, levels = c(tmp$feature))
  
  modified_y_labels <- sub("_[^_]*$", "", tmp$feature)
  
  # plot
  plot = tmp %>%
    ggplot(aes(x = feature, y = value, color = view)) +
    geom_hline(yintercept = 0, color = "grey", linetype = "dashed") +
    geom_point(size = 1.5) +  
    geom_segment(aes(xend = feature, yend = 0), linetype = "solid") +
    coord_flip() +  
    labs(title = factor_select, x = "", y = "Feature weight", color = "") +  
    scale_color_carto_d(palette = "Earth") +
    theme_bw() +
    scale_x_discrete(labels = setNames(modified_y_labels, tmp$feature))
  
  return(plot)
  
}
  
  