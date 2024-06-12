#' @name pcaBiplot
#' @description
#' This function creates a PCA biplot from a FactoMineR PCA object and corresponding data frame.
#' Customisable to plot varuing numbers of contributors, colour by different characteristics in supplied dataframe,
#' and whether to show individual lines.
#' @param pc.object a res.pca object output from FactoMineR
#' @param df.wide Wide dataframe where column 'primary' corresponds to rownames of the pc.object$ind$coord rownames.
#' @param colour.by Column to colour Biplot by, default is visitId
#' @param ncontrib Number of contributing variables to plot
#' @param individual.path Boolean to indicate whether individual plots should be shown or not, True by default
#' @param colours colours for colouring in the dots. Default values supplied and recommended.
#' @param dim.x Principal component to plot on axis x. PC1 by default
#' @param dim.y Principal component to plot on axis y. PCy by default
#' @return PCA biplot

pcaBiplot <- function(pc.object, df.wide,
                      colour.by = 'visitId',
                      ncontrib = 5,
                      individual.path = T,
                      colours = c('#48a0af', '#1b69a1', '#ec6669', '#f39668'),
                      dim.x = 1,
                      dim.y = 2){
  
  if(!identical(rownames(pc.object$ind$coord), df.wide$primary)){
    stop("PCA names and df names not identical; please check\n")
  }
  
  # Check whether both dimensions exist in PC object
  if(!all(c(paste0("Dim.", dim.x), paste0("Dim.", dim.y)) %in% colnames(pc.object$ind$coord))){
    stop("Not all dimensions found in PC object; please check\n")
  }
  
  # Generate 
  
  lab_pcx <- paste0('<b>PC', dim.x, '</b> (', round(pc.object$eig[dim.x, 2],2), '% of variance)')
  lab_pcy <- paste0('<b>PC2', dim.y, '</b> (', round(pc.object$eig[dim.y, 2],2), '% of variance)')
  
  plot <- factoextra::fviz_pca_biplot(pc.object,
                                      label = 'var',
                                      axes = c(dim.x, dim.y),
                                      addEllipses = F,
                                      pointshape = 19,
                                      col.var = 'grey20',
                                      labelsize = 3,
                                      col.ind = df.wide[[colour.by]],
                                      geom.ind="point", pointsize=1,
                                      select.var = list(contrib = ncontrib),
                                      repel = T)
  
  if(individual.path==T){
    plot <- plot +
      geom_path(aes(group = df.wide[['subjectId']]),
                colour = 'grey80',
                alpha = 0.9) +
      geom_point(aes(colour = df.wide[[colour.by]]),
                 size = 1,
                 alpha = 0.6) +
      stat_ellipse(aes(group = df.wide[[colour.by]],
                       colour = df.wide[[colour.by]])) +
      theme_bw() +
      guides(fill = element_blank()) +
      scale_colour_manual(name = 'month',
                          values = colours,
                          aesthetics = c('fill', 'colour')) +
      theme(axis.title.y = element_markdown(),
            axis.title.x = element_markdown(),
            legend.position = "inside",
            legend.position.inside = c(0.1, 0.2),
            legend.key.height= unit(0.5, 'cm'),
            legend.key.width= unit(0.5, 'cm'),
            legend.background = element_blank(),
            legend.text = element_text(size=10),
            legend.title=element_blank()) +
      labs(title = '',
           x = lab_pcx,
           y = lab_pcy)
  } else {
    plot <- plot +
      geom_point(aes(colour = df.wide[[colour.by]]),
                 size = 1,
                 alpha = 0.6) +
      stat_ellipse(aes(group = df.wide[[colour.by]],
                       colour = df.wide[[colour.by]])) +
      theme_bw() +
      guides(fill = element_blank()) +
      scale_colour_manual(name = 'month',
                          values = cols[c(2, 1, 4, 5)],
                          aesthetics = c('fill', 'colour')) +
      theme(axis.title.y = element_markdown(),
            axis.title.x = element_markdown(),
            legend.position = "inside",
            legend.position.inside = c(0.1, 0.2),
            legend.key.height= unit(0.5, 'cm'),
            legend.key.width= unit(0.5, 'cm'),
            legend.background = element_blank(),
            legend.text = element_text(size=10),
            legend.title=element_blank()) +
      labs(title = '',
           x = lab_pcx,
           y = lab_pcy)
  }
  
  return(plot)
}
