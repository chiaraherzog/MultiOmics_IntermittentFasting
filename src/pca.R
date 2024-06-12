#' @name pca
#' @description
#' A function to run a PCA and return data and plots
#' @param data A multiAssayExperiment object
#' @param experiment assay to be extracted (if multiple, concatenate them using c(''))
#' @param features features to be extracted alongside assay data (for biplot or PCA heatmap, e.g.)
#' @param complete Should only complete cases be considered? F by default.
#' @param complete.n Number of visits for complete cases. 4 by default, but can also be 2.
#' @param intervention Intervention to be extracted (string)
#' @param relabel relabeling for assay columns, an optionally supplied dataframe. column `x` of this dataframe needs to correspond to names in the multiassay experiment, with a column `label` with substituting values
#' @param ncp number of PCAs to include
#' @param ncontrib number of contributing features to include in the biplot, 5 by default
#' @return list object with data and plots.

pca <- function(data, experiment, features,
                complete = F, complete.n = 4,
                intervention = c('I', 'K'),
                relabel = NULL,
                npc = 5,
                ncontrib = 5){
  
  # Extract experiment data in wide format
  wide <- as.data.frame(wideFormat(data[,data@colData$interventionId %in% intervention & !data@colData$visitId %in% c("M12", "M18"),
                                        experiment],
                                   colData = features))
  
  # Wide format appends assay - we want to remove this:
  if(length(experiment) == 1){
    experiment_rowname <- paste0(gsub(":", ".", gsub(" ", ".", experiment)), "_")
    colnames(wide) <- gsub(experiment_rowname, "", colnames(wide))
  } else {
    for (i in experiment){
      experiment_rowname <- paste0(gsub(":", ".", gsub(" ", ".", i)), "_")
      colnames(wide) <- gsub(experiment_rowname, "", colnames(wide))
    }
  }
  
  
  # Optional: relabel variables
  if(!is.null(relabel)){
    wide <- wide |> 
      dplyr::select(any_of(c( 'primary', features, relabel$x))) |> 
      dplyr::rename_at(vars(relabel[relabel$x %in% colnames(wide),]$x), ~relabel[relabel$x %in% colnames(wide),]$label) 
  }
  
  # Optional: filter complete cases only
  if(complete == T){
    completecases <- wide |>
      dplyr::group_by(subjectId) |>
      dplyr::count() |>
      dplyr::filter(n == complete.n)
    
    wide <- wide |> 
      dplyr::filter(subjectId %in% completecases$subjectId)
    rm(completecases)
  }
  
  # Extract data only (i.e. removing covariates/features) for PCA
  data_pca <- wide |> 
    dplyr::select(-all_of(features)) |> 
    tibble::column_to_rownames('primary') |> 
    dplyr::select_if(~ !any(is.na(.)))
  
  # Run PCA
  res.pca <- FactoMineR::PCA(data_pca, ncp = 10, graph = F)
  
  # A: Scree diagram -----------------
  scree <- res.pca$eig
  
  scree <- scree |> 
    as.data.frame() |> 
    dplyr::slice(1:npc) |> 
    ggplot(aes(x = 1:npc,
               y = `percentage of variance`)) +
    geom_col(fill = '#1b69a1') +
    theme_bw() +
    theme(axis.title.y = element_markdown(),
          panel.grid = element_blank(),
          axis.ticks.x = element_blank()) +
    labs(x= 'PC', y = 'Percentage of variance explained (%)') +
    scale_y_continuous(expand = c(0.1, 0.1)) +
    scale_x_continuous(breaks = c(1:5),
                       labels = c(1:5))
  
  # B: Biplot ----------------
  biplot <- pcaBiplot(res.pca, wide,
                      ncontrib = ncontrib)
  
  # C: PCA Heatmap ----------------
  heatmap <- pcaHeatmap(res.pca, wide, features = features)
  
  
  return(list(pc.object = res.pca,
              wide = wide,
              scree = scree,
              biplot = biplot,
              heatmap = heatmap))
}
