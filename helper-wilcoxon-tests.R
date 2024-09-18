library(dplyr)
library(broom)
library(here)
library(MultiAssayExperiment)
here::i_am("helper-wilcoxon-tests.R")

load("data/data_raw.Rdata")
load(here("src/vars.Rdata"))
vars <- vars |>
  dplyr::filter(!grepl("ASVs$|families_clr", assay))


# extract all data in long format plus basic anno
df <- as.data.frame(longFormat(data, colData = c("age_at_consent", "interventionId", "subjectId", "visitId", "time", 'compliance', 'comprate', 'supplement',
                                                 'bmi_at_consent',
                                                 'pattern',
                                                 'mpstatrs'))) |> 
  dplyr::filter(!is.na(value)) |> 
  dplyr::inner_join(vars, by = c("rowname" = "x",
                                 "assay" = "assay")) |> 
  dplyr::select(-c(label, assay2)) |> 
  dplyr::mutate(rowname = paste0(assay, "_", rowname))
# all(df$rowname %in% vars$x)

# remove S and any non M0-M6
df_raw <- df |>
  dplyr::filter(interventionId!='S' & !visitId %in% c("M12", "M18")) |>
  dplyr::mutate(compliance = factor(compliance, levels = c("low", 'medium', 'high')),
                mpstatrs = factor(mpstatrs, levels = c("no", "yes"))) 

# remove full data
rm(data);gc()

# load in data
load("data/data_baseline_change.Rdata")

# extract all data in long format plus basic anno
df <- as.data.frame(longFormat(data, colData = c("age_at_consent", "interventionId", "subjectId", "visitId", "time", 'compliance', 'comprate', 'supplement',
                                                 'bmi_at_consent', 'pattern', 'mpstatrs'))) |> 
  dplyr::filter(!is.na(value)) |> 
  dplyr::inner_join(vars, by = c("rowname" = "x",
                                 "assay" = "assay")) |> 
  dplyr::select(-c(label, assay2)) |> 
  dplyr::mutate(rowname = paste0(assay, "_", rowname))

# remove S and any non M0-M6
df_change <- df |>
  dplyr::filter(interventionId!='S' & !visitId %in% c("M12", "M18")) |>
  dplyr::mutate(compliance = factor(compliance, levels = c("low", 'medium', 'high')),
                mpstatrs = factor(mpstatrs, levels = c("no", "yes")))

# remove full data
rm(data, df);gc()

source(here("src/wilcoxon_tests.R"))
wilcoxon_tests <- wilcoxon_tests(df_change, df_raw)

save(wilcoxon_tests, file = here("out/wilcoxon-tests.Rdata"))