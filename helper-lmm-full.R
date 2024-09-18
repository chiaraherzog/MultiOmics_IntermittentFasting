## ----setup, include=F-------------------------------------------------------------------------------------------------------------------------------------------------
## knitr::opts_chunk$set(echo = T, message = F, warning = F)


## ----libs, eval = T---------------------------------------------------------------------------------------------------------------------------------------------------
library(dplyr)
library(lmerTest)
library(lme4)
library(here)
library(broom)
library(MultiAssayExperiment)

here::i_am("helper-lmm-full.R")
source("src/summary_lmm_v2.R")


## ----getdata, eval = T------------------------------------------------------------------------------------------------------------------------------------------------
# load in data
load("data/data_normalized.Rdata")

# extract all data in long format plus basic anno
df <- as.data.frame(longFormat(data, colData = c("age_at_consent", "interventionId", "subjectId", "visitId", "time",
                                                 'compliance', 'comprate', 'supplement',
                                                 'bmi_at_consent',
                                                 'dbmi',
                                                 'mpstatrs'))) |> 
  dplyr::filter(!is.na(value))


# remove S and any non M0-M6
df <- df |>
  dplyr::filter(interventionId!='S' & !visitId %in% c("M12", "M18")) |>
  dplyr::mutate(compliance = factor(compliance, levels = c("low", 'medium', 'high')),
                visitId = factor(visitId, levels = c("M0", "M2", "M4", "M6")),
                visitIdOrdered = ordered(visitId)) |> 
  dplyr::mutate(rowname = paste0(assay, "_", rowname))

# remove full data
rm(data);gc()

data  <- df
variables = unique(data$rowname)
variablesko <- variables[grepl(": KO", variables)]


cat('Starting model analysis...\n')
## ----eval = T---------------------------------------------------------------------------------------------------------------------------------------------------------
# out_lmm <- summary_lmm(data,
#                        variables = variables,
#                        timeType = 'factor',
#                        outName = 'out_lmm_factor.Rdata')

out_lmm <- summary_lmm(data,
                       variables = variablesko,
                       timeType = 'factor',
                       outName = 'out_lmm_factor_ko.Rdata')


## ----lmm.cont, eval = F-----------------------------------------------------------------------------------------------------------------------------------------------
## out_lmm <- summary_lmm(data,
##                        variables = variables,
##                        timeType = 'continuous',
##                        outName = 'out_lmm_continuous_15Apr2024.Rdata')


## ----save.as.script, eval = F-----------------------------------------------------------------------------------------------------------------------------------------
## knitr::purl("helper-lmm.qmd",
##             output = 'helper-lmm.R')


## ----eval = F---------------------------------------------------------------------------------------------------------------------------------------------------------
## load("out/out_lmm_factor_15Apr2024.Rdata")
## 
## vars_clin <- read.table("src/clinical_variables.csv", header = T, sep = ',')
## 
## vars_cerv_meth <- readxl::read_xlsx("src/indices.xlsx", sheet = 1) |>
##   dplyr::add_row(x = 'ic',
##                  label = 'immune cell proportion')
## 
## vars_buccal_meth <- readxl::read_xlsx("src/indices.xlsx", sheet = 2) |>
##   dplyr::add_row(x = 'ic',
##                  label = 'immune cell proportion')
## 
## vars_blood_meth <- readxl::read_xlsx("src/indices.xlsx", sheet = 3)
## 
## load("src/populations_names_annotated.Rdata")
## vars_imm <- populations |>
##   dplyr::filter(main_analysis == 'yes') |>
##   dplyr::mutate(x = name,
##                 label = `population name`)
## 
## vars_lut <- plyr::rbind.fill(vars_clin, vars_cerv_meth, vars_buccal_meth, vars_blood_meth, vars_imm)
## 
## # Filter and rename vars
## n <- names(out_lmm)
## lmm2 <- lapply(n, function(x){
## 
##   out_lmm[[x]] |>
##     tidyr::separate('x', "_",
##                     into = c("assay", "x"),
##                     extra = "merge") |>
##     dplyr::left_join(vars_lut) |>
##     dplyr::mutate(variable = ifelse(!is.na(label), label, x)) |>
## 
##     # Filter superfluous assays
##     dplyr::filter(!grepl("Immune age: IF group|Immune age: SMK group|resonance$|stimulated", assay) & !grepl("ImmAge_gen_adj", assay)) |>
## 
##     # Filter methylation + clinical features with NAs in variabke
##     dplyr::filter((grepl("methylation|haemogram|composition|sonography|loss|cytometry", assay) & !is.na(label)) | !grepl("methylation|haemogram|composition", assay)) |>
## 
##     dplyr::select(-starts_with("std.error_")) |>
##     dplyr::arrange(assay) |>
##     dplyr::relocate(assay, variable) |>
##     dplyr::distinct()
## 
## 
## })
## names(lmm2) <- n
## 
## 
## 
## 
## n <- names(t_test)
## 
## t_test <- lapply(n, function(x){
##   t_test[[x]] |>
##     dplyr::left_join(vars_lut, by = c('variable' = 'x')) |>
##     dplyr::mutate(variable = ifelse(!is.na(label), label, variable)) |>
##     dplyr::select(-any_of(c("x", "label"))) |>
##     dplyr::arrange(assay, pick(any_of(c("p_M6", "p value (high versus other)", "p value (I versus K, high compliance)"))))
## })
## names(t_test) <- n
## 
## 
## # Rename additional assays to be more informative + tidy up formating and columns
## t_test <- lapply(n, function(x){
##   t_test[[x]] |>
##     dplyr::mutate(assay = ifelse(!is.na(assay2), assay2, assay)) |>
##     dplyr::select(-any_of(c("population name", "name", "fixable viability dye and antibodies α-", "staining", "type", "antibodies", "targets", "main_analysis",
##                             "assay2"))) |>
##     dplyr::mutate(across(contains(c("p_M2", "p_M4", "p_M6", "p value (high versus other)", "p value (I versus K, high compliance)")), ~ signif(., digits = 3)))
## })
## names(t_test) <- n
## 


## ----eval = F---------------------------------------------------------------------------------------------------------------------------------------------------------
## DT::datatable(lmm2[["Minimal model"]])


## ----eval = F---------------------------------------------------------------------------------------------------------------------------------------------------------
## DT::datatable(lmm2[["Basic model with BMI"]])


## ----eval = F---------------------------------------------------------------------------------------------------------------------------------------------------------
## DT::datatable(lmm2[["Intervention (higher compliance)"]])


## ----writelmm, eval = F-----------------------------------------------------------------------------------------------------------------------------------------------
## for (x in c("Minimal model", "Basic model with BMI", "Intervention (higher compliance)")){
##   file = case_when(x == 'Minimal model' ~ "6",
##                    x == 'Basic model with BMI' ~ "7",
##                    TRUE ~ "8")
## 
##   writexl::write_xlsx(as.data.frame(lmm2[[x]]),
##                       path = paste0("out/Extended-Data-Table-",
##                                     file,
##                                     ".xlsx")
##     )
## }

