## helper functions for bottum-up FDR testing lmm results microbiome
# Author: Charlotte Vavourakis

# use prepDataBouth and calcBouth for ASV-level results
# use prepDataBouth2 and calcBouth2 for family-level results
# use filterSignBouth to prep data in long format only with features that were significant after bottom up FDR correction

# workflow depends on makeCustomASVtab to make an ASV table from wide data with phylum, class, order, family tax info added

prepDataBouth <- function(datLmmOut, exp, datASV, estimate, prevalence){
  
  require(tidyverse)
  
  out = c()
  
  datLmmOut = datLmmOut %>% filter(str_detect(x, exp))
  datLmmOut$x = gsub(paste0(exp,"_"),"",datLmmOut$x)
  
  ## vector with p-values
  extract <- paste0("p.value_",estimate)
  pval.leaves <- datLmmOut %>% pull(extract)
  names(pval.leaves) = datLmmOut %>% pull(x)
  
  ## tax table
  taxTab = datASV %>% 
    select(Kingdom,Phylum,Class,Order,Family,OTU) %>%
    dplyr::rename(ASV=OTU) %>%
    filter(ASV %in% names(pval.leaves)) %>%
    arrange(Family, Order, Class, Phylum, Kingdom)
  taxTab$Kingdom = paste0("k__",taxTab$Kingdom)
  taxTab$Phylum = paste0("p__",taxTab$Phylum)
  taxTab$Class = paste0("c__",taxTab$Class)
  taxTab$Order = paste0("o__",taxTab$Order)
  taxTab$Family = paste0("f__",taxTab$Family)
  
  ## keep only prevalent ASVs
  datASV = datASV %>%
    filter(OTU %in% c(taxTab$ASV)) %>%
    select(-Kingdom,-Phylum,-Class,-Order,-Family,-Genus,-Species) %>%
    column_to_rownames(var="OTU")
  non_zero_counts = apply(datASV, 1, function(row) sum(row != 0))
  prevalence_df = data.frame(Feature = rownames(datASV),
                             Prevalence = (non_zero_counts / ncol(datASV)) * 100)
  keep = prevalence_df %>%
    filter(!Prevalence < prevalence) %>%
    pull(Feature)
  
  taxTab = taxTab %>% filter(ASV %in% keep)
  pval.leaves <- pval.leaves[keep]
  
  ## check order
  order = match(taxTab$ASV, names(pval.leaves))
  pval.leaves = pval.leaves[order]
  
  ## make list
  out[[1]] = taxTab
  out[[2]] = pval.leaves
  names(out) = c("tax.tab","pval.leaves")
  
  return(out)
}

calcBouth <- function(datLmmOut, exp, datASV, estimates, prevalence) {
  
  require(tidyverse)
  require(BOUTH)
  
  resultsBouth = list()
  
  for (i in 1:length(estimates)){
    datIn = prepDataBouth(datLmmOut, exp, datASV, estimates[i], prevalence)
    tryCatch({
      test = bouth(anno.table = datIn[[1]], pvalue.leaves = datIn[[2]],
                   na.symbol = "unknown", far = 0.1, is.weighted = TRUE)
      
      resultsBouth[[i]] = test$results.by.node %>%
        filter(is.detected == TRUE) %>%
        mutate(estimate = estimates[i])
    }, error = function(e) {
      cat("Error occurred for estimate:", estimates[i], "\n")
      resultsBouth[[i]] = NULL
    })
  }
  
  resultsBouth = do.call(rbind,resultsBouth)
  
  return(resultsBouth)
  
}

prepDataBouth2 <- function(datLmmOut, exp, datASV, estimate){
  
  require(tidyverse)
  
  out = c()
  
  datLmmOut = datLmmOut %>% filter(str_detect(x, exp))
  datLmmOut$x = gsub(paste0(exp,"_"),"",datLmmOut$x)
  
  ## vector with p-values
  extract <- paste0("p.value_",estimate)
  pval.leaves <- datLmmOut %>% pull(extract)
  names(pval.leaves) = datLmmOut %>% pull(x)
  
  ## tax table
  taxTab = datASV %>% 
    select(Kingdom,Phylum,Class,Order,Family) %>%
    distinct() %>%
    filter(Family %in% names(pval.leaves)) %>%
    arrange(Order, Class, Phylum, Kingdom)
  taxTab$Kingdom = paste0("k__",taxTab$Kingdom)
  taxTab$Phylum = paste0("p__",taxTab$Phylum)
  taxTab$Class = paste0("c__",taxTab$Class)
  taxTab$Order = paste0("o__",taxTab$Order)
  
  ## check order
  order = match(taxTab$Family, names(pval.leaves))
  pval.leaves = pval.leaves[order]
  
  ## make list
  out[[1]] = taxTab
  out[[2]] = pval.leaves
  names(out) = c("tax.tab","pval.leaves")
  
  return(out)
}

calcBouth2 <- function(datLmmOut, exp, datASV, estimates) {
  
  require(tidyverse)
  require(BOUTH)
  
  resultsBouth = list()
  
  for (i in 1:length(estimates)){
    datIn = prepDataBouth2(datLmmOut, exp, datASV, estimates[i])
    tryCatch({
      test = bouth(anno.table = datIn[[1]], pvalue.leaves = datIn[[2]],
                   na.symbol = "unknown", far = 0.1, is.weighted = TRUE)
      
      resultsBouth[[i]] = test$results.by.node %>%
        filter(is.detected == TRUE) %>%
        mutate(estimate = estimates[i])
    }, error = function(e) {
      cat("Error occurred for estimate:", estimates[i], "\n")
      resultsBouth[[i]] = NULL
    })
  }
  
  resultsBouth = do.call(rbind,resultsBouth)
  
  return(resultsBouth)
  
}

filterSignBouth <- function(resultsBouth, estimatesPlot,datASV_custom){
  
  # filter Bouth results for what to plot
  sign = resultsBouth %>% 
    filter(estimate %in% estimatesPlot) %>%
    mutate(taxget = case_when(
      ASV != " " ~ ASV,
      ASV == " " & Family != " " ~ Family,
      ASV == " " & Family == " " & Order != " " ~ Order,
      ASV == " " & Family == " " & Order == " " & Class != " " ~ Class,
      ASV == " " & Family == " " & Order == " " & Class == " " & Phylum != " " ~ Phylum
    )) %>%
    mutate(taxget2 = gsub(".*?__", "", taxget)) %>%
    na.omit() # do not consider Bacteria overall
  
  taxget = unique(sign$taxget)
  
  abunddat = list()
  
  # get changes in relative abundance from baseline, sum for taxa above ASV-level
  for (i in 1:length(taxget)){
    if(grepl("ASV", taxget[i])) {
      
      abunddat[[i]] = datASV_custom %>% 
        filter(ASV == taxget[i]) %>%
        select(-Phylum,-Class,-Order,-Family) %>%
        dplyr::rename(taxget2=ASV)
      
    } else if (grepl("f__", taxget[i])) {
      
      abunddat[[i]] = datASV_custom %>% 
        filter(Family == gsub(".*?__", "", taxget[i])) %>%
        group_by(Family) %>%
        summarize(across(contains("SM"), sum, na.rm = TRUE)) %>%
        dplyr::rename(taxget2=Family) %>%
        as.data.frame()
      
    } else if (grepl("o__", taxget[i])) {
      
      abunddat[[i]] = datASV_custom %>% 
        filter(Order == gsub(".*?__", "", taxget[i])) %>%
        group_by(Order) %>%
        summarize(across(contains("SM"), sum, na.rm = TRUE)) %>%
        dplyr::rename(taxget2=Order) %>%
        as.data.frame()
      
    } else if (grepl("c__", taxget[i])) {
      
      abunddat[[i]] = datASV_custom %>% 
        filter(Class == gsub(".*?__", "", taxget[i])) %>%
        group_by(Class) %>%
        summarize(across(contains("SM"), sum, na.rm = TRUE)) %>%
        dplyr::rename(taxget2=Class) %>%
        as.data.frame()
      
    } else if (grepl("p__", taxget[i])) {
      
      abunddat[[i]] = datASV_custom %>% 
        filter(Phylum == gsub(".*?__", "", taxget[i])) %>%
        group_by(Phylum) %>%
        summarize(across(contains("SM"), sum, na.rm = TRUE)) %>%
        dplyr::rename(taxget2=Phylum) %>%
        as.data.frame()
      
    }
  }
  
  abunddat = do.call(rbind,abunddat)
  
  dat = full_join(sign,abunddat)
  
  return(dat)
}

