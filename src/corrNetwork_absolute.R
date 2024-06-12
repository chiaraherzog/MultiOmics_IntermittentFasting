corrNetwork_absolute <- function(data,
                        pathwayItems,
                        experiments,
                        corMethod = 'pearson',
                        minCor = 0.2,
                        visit){
  
  tmp <- longFormat(data[, data$visitId==visit & data$interventionId != 'S', experiments],
                    colDataCols = "subjectId") |>
    as.data.frame() |> 
    dplyr::mutate(label = paste(assay, rowname, sep = "_")) |> 
    dplyr::filter(grepl(pathwayItems, label, ignore.case = T)) |> 
    tidyr::pivot_wider(id_cols = c("subjectId"),
                       values_from = 'value',
                       names_from = label) |> 
    tibble::column_to_rownames('subjectId')
  
  vars_relabel <- data.frame(name = colnames(tmp)) |> 
    tidyr::separate(name, "_", into = c("group", "label"), extra = 'merge', remove = F)
  
  g <- tmp |> 
    corrr::correlate(method = corMethod) |> 
    corrr::stretch() |> 
    dplyr::filter(abs(r) > minCor) |> 
    tidygraph::as_tbl_graph(directed = FALSE) |> 
    tidygraph::activate(nodes) |> 
    dplyr::left_join(vars_relabel) |>
    tidygraph::activate(edges) |> 
    dplyr::mutate(weight=r,
                  assoc = ifelse(r < 0, 'neg', 'pos'))
  
  set.seed(1)
  
  p <- g |> 
    ggraph(layout = "graphopt") + 
    geom_edge_link(aes(width = abs(weight),
                       color = assoc), alpha = 0.2) + 
    scale_edge_width(range = c(0.2, 1.5),
                     name = 'absolute r') +
    scale_colour_manual(values = cols,
                        name = '') +
    scale_edge_color_manual(values = c("gray70", cols[5]),
                            name = 'sign',
                            labels = c('negative',
                                       'positive')) +
    geom_node_point(aes(colour = group),
                    size = 4) +
    geom_node_text(aes(label = label), repel = TRUE,
                   size = 2.7) +
    theme_graph() +
    theme(legend.title = element_markdown())
  
  print(p)
  
}
