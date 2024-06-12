network_diag <- function(measure){
  
  library(GGally)
  library(network)
  library(sna)
  library(ggplot2)
  library(igraph)
  
  inner <- corr |> 
    dplyr::filter(if_any(c('measure1', 'measure2'), ~ grepl(measure, .)))
  
  vars <- unique(c(inner$measure1, inner$measure2))
  load(here("src/populations_names_annotated.Rdata"))
  
  # now let's expand the network:
  expanded <- corr |> 
    dplyr::filter(if_any(c('measure1', 'measure2'), ~ grepl(paste0(vars, collapse = "|"), .))) |> 
    tidyr::separate(measure1, "_", into = c(NA, 'measure1'), extra = 'merge') |> 
    tidyr::separate(measure2, "_", into = c(NA, 'measure2'), extra = 'merge') |> 
    dplyr::left_join(dplyr::select(populations, name, `population name`), by = c('measure1' = 'name')) |> 
    dplyr::left_join(dplyr::select(populations, name, `population name`), by = c('measure2' = 'name')) |> 
    dplyr::mutate(measure1 = ifelse(is.na(`population name.x`), measure1, `population name.x`)) |> 
    dplyr::mutate(measure2 = ifelse(is.na(`population name.y`), measure2, `population name.y`))
  
  g <- graph_from_edgelist(as.matrix(expanded[,1:2]), directed = F)
  vertices <- names(V(g))
  assays1 <- expanded[match(vertices, c(expanded$measure1)),]$assay1
  assays2 <- expanded[match(vertices, c(expanded$measure2)),]$assay2
  assays1[is.na(assays1)] <- assays2[is.na(assays1)]
  assays1 <- gsub("-log", "", assays1)
  cols <- grid.col[match(assays1, names(grid.col))]
  
  V(g)$color <- cols
  p <- plot(g, vertex.color=V(g)$color,
       vertex.label.color = 'black', vertex.label.family = 'Helvetica',
       vertex.size = 10,
       vertex.label.cex = 0.7,
       layout = layout.fruchterman.reingold)
  
  return(p)
}