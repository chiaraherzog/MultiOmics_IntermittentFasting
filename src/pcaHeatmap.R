#' @name pcaHeatmap
#' @description
#' This function creates a PCA heatmap from a FactoMineR PCA object and corresponding data frame.
#' @param pc.object a res.pca object output from FactoMineR
#' @param df.wide Wide dataframe where column 'primary' corresponds to rownames of the pc.object$ind$coord rownames.
#' @param npc Number of principal components to include in the heatmap
#' @param features to include
#' @return PCA heatmap

pcaHeatmap <- function(pc.object, df.wide,
                       features = features,
                       npc = 5){
  
  if(!identical(rownames(pc.object$ind$coord), df.wide$primary)){
    stop("PCA names and df names not identical; please check\n")
  }
  
  # Select PCs
  pcs <- pc.object$ind$coord |>
    as.data.frame() |> 
    dplyr::select(Dim.1:paste0("Dim.", npc))
  
  # Select Pheno; remove visitId as duplicated (time should be included)
  pheno <- df.wide |> 
    dplyr::select(any_of(features)) |> 
    dplyr::mutate_if(is.character, as.factor)
  
  # If mpstatrs, remove stars for simplification
  if('mpstatrs' %in% colnames(pheno)){
    pheno$mpstatrs <- gsub("[*]", "", pheno$mpstatrs)
  }
  
  # Compute associations
  mat <- matrix(ncol = ncol(pcs),
                nrow = ncol(pheno))
  colnames(mat) <- colnames(pcs)
  rownames(mat) <- colnames(pheno)
  
  for (i in colnames(mat)){
    for (j in rownames(mat)){
      
      if(is.numeric(pheno[,j])){
        mat[j,i] <- cor.test(pcs[,i], pheno[,j])$p.value
      } else {
        mat[j,i] <- kruskal.test(pcs[,i], pheno[,j])$p.value
      }
      
    }
  }
  
  map <- mat
  map <- apply(map, 2, function(t) ifelse(t < 0.05, t, NA))
  
  # Relabel columns and rows for better interpretability
  labLookup =  c('interventionId' = 'interventionId',
                 'subjectId' = 'subjectId',
                 'visitId' = 'visitId',
                 'time' = 'visitId (numerical)',
                 'compliance' = 'visitId (categorical)',
                 'comprate' = 'compliance rate (%)',
                 'smoking_ever' = 'smoking (ever)',
                 'mpstatrs' = 'menopausal status (at consent)',
                 'age_at_consent' = 'age (at consent)',
                 'bmi_at_consent' = 'BMI (at consent)', 
                 'etohu_curr' = 'current alcohol units/wk',
                 'diet' = 'dietary preference (at consent)',
                 'pattern' = 'fasting pattern',
                 'preg_ever' = 'pregnancy (ever)',
                 'ocp_curr' = 'OCP (at consent)',
                 'hrt_curr' = 'HRT (at consent)',
                 'intactcurr' = 'intense activity (current/week)')
  
  ind <- match(rownames(map), names(labLookup))
  rownames(map) <- labLookup[ind]
  
  colnames(map) <- gsub("Dim.", "PC", colnames(map))
  
  # plot the heatmap
  pcaheatmap <- Heatmap(-log10(map[,1:npc]),
                        cluster_columns = F,
                        row_names_side = 'left',
                        cluster_rows = F,
                        column_title = NULL,
                        show_row_dend = F,
                        show_column_dend = F,
                        name = '-log10(p value)',
                        na_col = 'white',
                        col = circlize::colorRamp2(breaks = seq(40, 1.3, length.out = 5),
                                                   colors = rev(viridis::viridis(5)),
                        ),
                        row_names_gp = grid::gpar(fontsize = 9),
                        column_names_gp = grid::gpar(fontsize = 9),
                        border_gp = gpar(lwd = 0.5),
                        border = T)
  return(pcaheatmap)
}
