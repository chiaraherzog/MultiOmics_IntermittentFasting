# plot relative abundances from a ASVtable (counts)
# Author: Charlotte Vavourakis
# version for IF

plot_taxa_heatmap <- function(ASVtable, pheno, agg, add, abu, type, colors){
  
  # to do fix colors?
  
  require(ampvis2)
  require(tidyverse)
  require(ComplexHeatmap)
  require(rcartocolor)
  
  # theme_set(theme_minimal(base_size=8))
  # 
  ## order pheno2 for plots without hierarchical clustering
  
  # load data for analysis with ampvis2, convert to rel abundance
  d = amp_load(otutable=ASVtable, metadata=pheno)
  d = normaliseTo100(d)
  
  # aggregate
  d = aggregate_abund(
    d$abund,
    d$tax,
    tax_aggregate = agg,
    tax_add = add,
    format = "abund",
    calcSums = FALSE)
  
  # abundance filter
  d = d %>% 
    mutate(threshold_bool = apply((select(., ends_with(type))),1,function(x) +(any(x > abu)))) %>%
    filter(threshold_bool == 1) %>%
    select(-threshold_bool) %>%
    droplevels()

  # annotation
  vars = c("compliance",
           #"pattern",
           "interventionId",
           "visitId")
  
  annot = pheno %>%
    select(all_of(c("sampleId",vars))) %>%
    column_to_rownames(var="sampleId")
  
  annot$compliance <- factor(annot$compliance, levels = c("low","medium","high"))
  
  mycolors1 = cols[c(2, 1, 4, 5)]
  names(mycolors1) <- unique(annot$visitId)
  
  mycolors2 = cols[c(2,7)]
  names(mycolors2) <- unique(annot$interventionId)
  
  mycolors3 = cols[c(4, 5, 1)]
  names(mycolors3) <- c("low","medium","high")
  
  # mycolors4 = viridis_pal(option = "turbo")(5)
  # names(mycolors4) <- c(unique(annot$pattern),"NA")
  
  mycolors <- list(interventionId = mycolors2,
                   compliance = mycolors3,
                   #pattern = mycolors4,
                   visitId = mycolors1)
  
  # heatmap
  plot = pheatmap(d, 
                  annotation_col=annot,
                  annotation_colors = mycolors,
                  cluster_cols = FALSE,
                  show_colnames=FALSE,
                  border_color = NA)
  
  return(plot)
}