
## PCA with biplot and covariate association plots for the microbiome data
# Hellinger transformation and scaling = F
# version for IF
# Author: Charlotte Vavourakis

pc_microbiomes <- function(df_raw_wide,features,vegan_trans, n_pc, sampletype, label_col){
  
  require(tidyverse)
  require(vegan)
  require(factoextra)
  require(FactoMineR)
  require(ComplexHeatmap)
  require(viridis)
  
  plotList <- list()
  
  ### perform PCA ###----
  n1 <- ncol(df_raw_wide) - length(features)
  n2 <- n1 + 1
  
  pcmat <- df_raw_wide[,c(1:n1)] |> 
    tibble::column_to_rownames('primary') |> 
    select_if(~ !any(is.na(.)))
  
  if (vegan_trans == "hellinger") {
    pcmat <- vegan::decostand(pcmat, method = "hellinger")
    res.pca <- FactoMineR::PCA(pcmat, scale=F, ncp=10, graph = FALSE) 
  } else {
    res.pca <- FactoMineR::PCA(pcmat, scale=T, ncp=10, graph = FALSE) 
  }
  
  ### Screeplot ###----
  plotList[[1]] <- factoextra::fviz_screeplot(res.pca,
                                              addlabels = F,
                                              linecolor = NA,
                                              barcolor = cols[1],barfill = cols[1]) +
    theme_bw() +
    labs(x = "PC", y = "Percentage of variance explained (%)", subtitle = "", title = "") +
    theme(aspect.ratio = 1.5)
  
  lab_pc1 <- paste0('PC1 (', round(res.pca$eig[1,2]), '% of variance)')
  lab_pc2 <- paste0('PC2 (', round(res.pca$eig[2,2]), '% of variance)')
  
  
  ### PC1 vs PC2 ###----
  plotList[[2]] <- factoextra::fviz_pca_biplot(res.pca,
                                               label = 'var',
                                               axes = c(1, 2),
                                               habillage = as.factor(df_raw_wide$visitId),
                                               addEllipses = T,
                                               geom = c("point", "arrow", "text"),
                                               pointshape = 19,
                                               geom.ind="point", pointsize=1, alpha.ind = 0.5,
                                               select.var = list(contrib = 3),
                                               repel = T,
                                               ellipse.alpha = 0,
                                               col.var = "black") +
    theme_bw() +
    theme(aspect.ratio = 1) +
    guides(fill = element_blank()) +
    scale_colour_manual(name = 'month',
                        values = cols[c(2, 1, 4, 5)],
                        aesthetics = c('fill', 'colour')) +
    labs(x = lab_pc1,y = lab_pc2) 
  
  ### Factor analysis  ###----
  pc_dat <- cbind(res.pca$ind$coord, df_raw_wide[df_raw_wide$primary %in% rownames(pcmat),c(1, n2:ncol(df_raw_wide))])
  colnames(pc_dat) <- gsub("Dim.","PC",colnames(pc_dat))
  
  pc_dat <- pc_dat |> 
    dplyr::filter(!is.na(compliance)) %>%
    dplyr::mutate(mpstatrs = gsub("[*]", "", mpstatrs)) 
  
  pcs <- pc_dat |> 
    dplyr::select(PC1:PC10)
  
  tmp <- pc_dat |> 
    dplyr::select(any_of(features)) |> 
    dplyr::mutate_if(is.character, as.factor)
  
  mat <- matrix(ncol = ncol(pcs),
                nrow = ncol(tmp))
  colnames(mat) <- colnames(pcs)
  rownames(mat) <- colnames(tmp)
  
  for (i in colnames(mat)){
    for (j in rownames(mat)){
      if(is.numeric(tmp[,j])){
        mat[j,i] <- cor.test(pcs[,i], tmp[,j])$p.value
        # pmat[j,i] <- cor.test(pcs[,i], tmp[,j])
      } else {
        mat[j,i] <- kruskal.test(pcs[,i], tmp[,j])$p.value
      }
    }
  }
  
  map <- mat
  map <- apply(map, 2, function(t) ifelse(t < 0.05, t, NA))
  map <- map[,1:n_pc]
  
  labs = c('interventionId', 'subjectId', 'visitId', 'compliance rate',
           'menopausal status', 'age (at consent)', 'BMI (at consent)', 
           'smoking (ever)', 'current alcohol units/wk', 'dietary pattern', 'pregnancy (ever)')
  
  plotList[[3]] <- Heatmap(-log10(map),
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
                           border = T,
                           top_annotation = HeatmapAnnotation(grp = anno_block(gp = gpar(fill = label_col,lwd = 0.5),
                                                                               labels = sampletype,
                                                                               labels_gp = gpar(col = "white",
                                                                                                fontsize = 10,
                                                                                                fontface = "bold")))
                           )
  
  return(plotList)
  
}
