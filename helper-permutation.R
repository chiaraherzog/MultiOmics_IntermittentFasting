## ----setup, include=F------------------------------------------------------------------------------------------------
knitr::opts_chunk$set(echo = T, message = F, warning = F, eval = F)


## ----libs, eval = T--------------------------------------------------------------------------------------------------
suppressMessages(suppressWarnings({
  library(dplyr)
  library(lmerTest)
  library(lme4)
  library(here)
  library(broom)
  library(MultiAssayExperiment)
}))

here::i_am("helper-permutation.qmd")
source("src/permutationTesting.R")


## ----getdata, eval = F-----------------------------------------------------------------------------------------------
# load in data
load("data/data_normalized.Rdata")

load("src/vars.Rdata")

varsList <- vars |> dplyr::filter(!assay %in% c('Saliva microbiome: ASVs', 'Stool microbiome: ASVs',
                                            'Saliva microbiome: families_clr', 'Stool microbiome: families_clr')) |>
  dplyr::group_by(assay) |>
  dplyr::summarise(features = list(x), .groups = 'drop') |> tibble::deframe()

# extract all data in long format plus basic anno
df <- as.data.frame(longFormat(subsetByAssay(data,varsList),
                             colData = c("age_at_consent", "interventionId", "subjectId", "visitId", "time",
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

# check class
class(df$visitId)
class(df$visitIdOrdered)

# remove full data
rm(data);gc()

data <- df
variables = unique(data$rowname)


## --------------------------------------------------------------------------------------------------------------------
cat("Starting basic")
outPerm <- permutationTesting(data, variables, n = 100)
save(outPerm, file = 'out/outPerm.Rdata')

cat("done.\nStarting minimal")

outPermMinimal <- permutationTesting(data, variables, n = 100,
                                     model = 'value ~ visitId + age_at_consent + bmi_at_consent + (1|subjectId)')
save(outPermMinimal, file = 'out/outPermMinimal.Rdata')

cat("done.\nStarting intervention")

# higher compliance
data <- data[data$compliance == 'high',]

outPermIntervent <- permutationTesting(data, variables, n = 100,
                                       model = 'value ~ visitId*interventionId + age_at_consent + bmi_at_consent + (1 | subjectId)')
save(outPermIntervent, file = 'out/outPermIntervention.Rdata')
cat("done.")

## ----save.as.script, eval = F----------------------------------------------------------------------------------------
# knitr::purl("helper-permutation.qmd",
#             output = 'helper-permutation.R')

