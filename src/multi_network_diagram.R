#' @name multi_network_diagram
#' @param features Features that should be included in the diagram
#' @param seed set for the layout
#' @param layout layout option
#' @param legend should a lgend be printed?

multi_network_diagram <- function(features,
                                  corrObj,
                                  seed = 7,
                                  layout = layout.fruchterman.reingold,
                                  legend = T, expand = T,
                                  expand_r = 0.3){
  
  # libraries
  library(GGally)
  library(network)
  library(sna)
  library(ggplot2)
  library(igraph)
  
  # Filter features relevant to the diagram
  inner <- corrObj |> 
    dplyr::filter(if_any(c('measure1', 'measure2'), ~ grepl(features, .)))
  
  if(expand == T){
    features_exp <- paste0(unique(c(inner$measure1, inner$measure2)), collapse = "|")
    inner <- corrObj |> 
      dplyr::filter(if_any(c('measure1', 'measure2'), ~ grepl(features, .) | 
                             (if_any(c('measure1', 'measure2'), ~ grepl(features_exp, .) & abs(r) >= expand_r))))
    features <- features_exp
  }
  
  features_upd <- unique(c(inner[grepl(features, inner$measure1),]$label1,
                           inner[grepl(features, inner$measure2),]$label2))
                           
  g <- graph_from_edgelist(as.matrix(inner[,c('label1', 'label2')]), directed = F)
  vertices <- names(V(g))
  assays1 <- inner[match(vertices, c(inner$label1)),]$assay1
  assays2 <- inner[match(vertices, c(inner$label2)),]$assay2
  assays1[is.na(assays1)] <- assays2[is.na(assays1)]
  assays1 <- gsub("-log", "", assays1)
  cols <- grid.col[match(assays1, names(grid.col))]
  
  col_edge <- ifelse(sign(inner$r) == -1, "lightblue", "coral3")
  names(col_edge) <- ifelse(sign(inner$r) == -1, "negative", "positive")
  rcorr <- scales::rescale(abs(inner$r), to = c(1, 8))
  names(rcorr) <- abs(inner$r)
  
  quantiles <- unname(quantile(rcorr, probs = c(0, 0.5, 0.8, 0.9, 1)))
  
  grid.col2 = grid.col[names(grid.col) %in% unique(gsub("-log", "", c(inner$assay1, inner$assay2)))]
    
  for (i in quantiles){
    
    if(i == quantiles[1]){
      x <- which.min(rcorr-i)
    } else {
      x <- c(x, which.min(abs(rcorr-i)))
    }
  }
  
  names(quantiles) <- signif(as.numeric(names(x)), 2)
  
  V(g)$color <- cols
  E(g)$color <- col_edge
  size_v = ifelse(names(V(g)) %in% features_upd, 15, 6)
  fontface_v = ifelse(names(V(g)) %in% features_upd, 2, 1)
  cex_v = ifelse(names(V(g)) %in% features_upd, 0.7, 0.6)
  
  # layout_nicely
  # layout_as_tree
  # layout_nicely
  
  set.seed(seed)
  p <- plot(g, vertex.color=V(g)$color,
            edge.color = E(g)$color,
            edge.width = rcorr,
       vertex.label.color = 'black',
       vertex.label.family = 'Helvetica',
       vertex.label.font = fontface_v,
       vertex.frame.color = NA,
       vertex.size = size_v,
       vertex.label.cex = cex_v,
       layout = layout,
       margin = 0)
  
  if(legend == T){
    legend("topleft",
           lty = 1,bty = "n",
           legend=unique(names(col_edge)),
           col=unique(col_edge), border=NA,
           title = "rmcorr")
    
    legend("left", bty = "n",
           lty = 1,
           legend=names(quantiles),
           lwd = quantiles,
           title = 'abs(correlation)')
    
    legend("bottomleft",bty = "n",
           legend=unique(gsub("\n", " ", names(grid.col2))),
           fill=unique(grid.col2), border=NA)
  }
  
  return(p)
}
