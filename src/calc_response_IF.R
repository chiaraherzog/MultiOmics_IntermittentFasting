
# Calculate excess weight loss in units of BMI (EBWL) at M2, M4, M6
# Calculate a score and score categories based on EBWL_M6 and comprate (compliance),
# reflecting how much weight is lost given the diet effort
# Author: Charlotte Vavourakis

calc_response_IF <- function(path_data = "~/Dropbox/eca/tirolgesund/7-manuscripts/if1/1-repository/data", 
                             samples_select=NULL){
  
  require(MultiAssayExperiment)
  require(dplyr)
  
  ## Calculate EBWL at different time points
  # gather bmi changes for IF study arm
  pathIn = file.path(path_data, "data_baseline_change.Rdata")
  tryCatch({
    load(pathIn)
  }, error = function(e) {
    stop(paste("Error: Unable to load file:", pathIn))
  })
  vars = "bmi"
  dat = as.data.frame(wideFormat(data[vars,,], 
                                 colData = c("interventionId","subjectId", "visitId", "bmi_at_consent","comprate"))) %>%
    filter(interventionId != "S") %>%
    select(-interventionId)
  
  # calculate floor based on absolute bmis recorded for IF study arm
  pathIn = file.path(path_data, "data_raw.Rdata")
  tryCatch({
    load(pathIn)
  }, error = function(e) {
    stop(paste("Error: Unable to load file:", pathIn))
  })
  vars = "bmi"
  floor = as.data.frame(wideFormat(data[vars,,], 
                                   colData = c("interventionId"))) %>%
    filter(interventionId != "S") %>%
    pull(Body.composition_bmi) %>%
    min()
  rm(data) ;gc()
  
  # EBWL calculation
  dat = dat %>%
    mutate(excessBMI_M0 = bmi_at_consent - floor) %>% 
    mutate(EBWL = Body.composition_bmi/excessBMI_M0 * 100 * (-1)) %>% # positive value means weight lost
    select(subjectId,visitId,EBWL,comprate)
  
  # rearrange to wide format, add means
  dat = dat %>%
    pivot_wider(names_from = visitId,
                names_glue = "{.value}_{visitId}",
                values_from = EBWL) %>%
    mutate(EBWL_mean = rowMeans(select(., starts_with("EBWL")))) %>%
    as.data.frame()
  
  # # sanity check with previous results
  # check <- readRDS("~/Dropbox/eca/tirolgesund/5-exploratory/5-microbiome/2-analysis-fasting_gut/0-output/pheno_stool_IK_full.Rds") %>%
  #   full_join(dat, by=c("subjectId")) %>%
  #   select(subjectId,EBWL_M6, M6_EBWL)
  
  ## calculate response score (residuals regression) based on EBWL_M6
  # continuous score
  dat2 = dat %>% select(subjectId,EBWL_M6,comprate) %>% distinct() %>% na.omit()
  model = lm(EBWL_M6 ~ comprate, data = dat2)
  dat2$predicted_weight_loss = predict(model, newdata = dat2) # predicted values
  dat2$response_IF = dat2$EBWL_M6 - dat2$predicted_weight_loss # calculate residuals and use as composite score
  
  # three equal categories across IF data for the score
  quantiles = quantile(dat2$response_IF, c(1/3, 2/3))
  dat2$response_IF_cat = cut(dat2$response_IF,
                             breaks = c(-Inf, quantiles[1], quantiles[2], Inf),
                             labels = c("Low Response", "Expected Response", "High Response"),
                             include.lowest = TRUE)
  
  dat2 = dat2 %>% select(subjectId,response_IF,response_IF_cat)
  dat = full_join(dat,dat2, by = "subjectId") %>% select(-EBWL_M0)
  
  # select subset of samples if user needs
  if (!is.null(samples_select)) {
    dat = dat %>% filter(subjectId %in% samples_select)
  }
  
  return(dat)
  
}