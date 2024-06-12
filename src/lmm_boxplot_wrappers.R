
# boxplot wrappers for plotting DA detected with lmm models
# depends on paired_longitudinal_compliance.R
# Author: Charlotte Vavourakis

boxplot_wrap_interaction <- function(dat1,dat2,exp,pval) {
  
  require(tidyverse)
  
  # select significant taxa
  toPlot = dat1 %>%
    dplyr::filter(`p.value_visitIdM6:interventionIdK` < pval | 
                    (`p.value_visitIdM6:interventionIdK` < 0.05 & `p.value_visitIdM6` < pval) ) %>%
    pull(x)
  
  # generate plots
  listPlots = list()
  for (i in 1:length(toPlot)){
    
    # default panel 1
    df = dat2 %>% filter(interventionId == "I")
    plota = paired_longitudinal_compliance(df,
                                           variable = toPlot[i],
                                           ylab = paste0(toPlot[i]), 
                                           p = 'p.signif') +
      ggtitle("intervention I")
    
    # default panel 2
    df = dat2 %>% filter(interventionId == "K")
    
    plotb = paired_longitudinal_compliance(df, 
                                           variable = toPlot[i],
                                           ylab = "", 
                                           p = 'p.signif') +
      ggtitle("intervention K")
    
    # Update plot y-axis limits
    y_min = min(plota$data$value, plotb$data$value)
    y_max = max(plota$data$value, plotb$data$value)
    y_max = round(y_max + 0.1 * y_max, digits=2)
    y_limits = c(y_min,y_max)
    plota = plota + coord_cartesian(ylim = y_limits)
    plotb = plotb + coord_cartesian(ylim = y_limits)
    
    listPlots[[i]] = plota + plotb
  }
  
  return(listPlots)
  
}

boxplot_wrap_overall <- function(dat1,dat2,exp,pval){
  
  require(tidyverse)
  
  # select significant taxa
  toPlot = dat1 %>%
    dplyr::filter(`p.value_visitIdM6` < pval) %>%
    dplyr::filter(!`p.value_visitIdM6:interventionIdK` < 0.05) %>%
    pull(x)
  
  # generate plots
  listPlots = list()
  for (i in 1:length(toPlot)){
    listPlots[[i]] = paired_longitudinal_compliance(dat2, 
                                                    variable = toPlot[i],
                                                    ylab = paste0(toPlot[i]), 
                                                    p = 'p.signif') # relative up
  }
  
  return(listPlots)
  
}

boxplot_wrap_fragmented_compliance <- function(dat1,dat2,exp,pval,visitId){
  
  # plot interaction term at specific time point for high compliance model
  
  require(tidyverse)
  
  filter_get <- paste0("p.value_visitId",visitId,":compliancehigh")
  
  # select significant taxa
  toPlot = dat1 %>%
    #dplyr::filter(`p.value_visitIdM6:compliancehigh` < pval) %>%
    dplyr::filter(!!sym(filter_get) < pval) %>%
    pull(x)
  
  # generate plots
  listPlots = list()
  for (i in 1:length(toPlot)){
    listPlots[[i]] = paired_longitudinal_compliance(dat2, 
                                                    variable = toPlot[i],
                                                    ylab = paste0(toPlot[i]), 
                                                    p = 'p.signif') # relative up
  }
  
  return(listPlots)
  
}
