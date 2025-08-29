# Which baseline features predict dbmi?


load("data/data_normalized.Rdata")

# extract all data in long format plus basic anno
df <- as.data.frame(longFormat(data, colData = c("age_at_consent", "interventionId", "subjectId", "visitId", "time",
                                                 'compliance', 'comprate', 'supplement',
                                                 'bmi_at_consent',
                                                 'dbmi',
                                                 'mpstatrs'))) |> 
  dplyr::filter(!is.na(value) & visitId == 'M0')

load("src/vars.Rdata")
# remove vars related to bmi
vars <- vars |> dplyr::filter(! x %in% c("weight", "fm"))

# remove S and any non M0-M6
df <- df |>
  dplyr::inner_join(vars, by = c('rowname' = 'x',
                                 'assay' = 'assay')) |> 
  dplyr::select(-c(label, assay2)) |> 
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

data  <- df
variables = unique(data$rowname)

median_dbmi <- median(data[data$visitId=='M0' & data$rowname=='Body composition_bmi', ]$dbmi, na.rm = T)
data$dbmigroup <- ifelse(data$dbmi <= median_dbmi, "higher weight loss", "lower weight loss")
data$dbmigroup <-factor(data$dbmigroup, levels = c("lower weight loss", "higher weight loss"))
data <- data[data$compliance=='high',]

cores = 3
library(parallel)
library(dplyr)
library(broom.mixed)
library(R.utils)

out <- mclapply(variables,  function(i) {
  tryCatch(lm(dbmi ~ bmi_at_consent + age_at_consent + interventionId + value ,
                data = data[data$rowname==i,]),
           error = function(e) return(e))
}, mc.cores = cores, mc.preschedule = T, mc.cleanup = T)
names(out) <- variables

## Remove any with error
out <- out |>
  purrr::keep( ~ !inherits(.x, 'error')) 

df <- dplyr::bind_rows(lapply(variables[variables %in% names(out)], function(x) as.data.frame(cbind(x, broom.mixed::tidy(out[[x]]))))) |> 
  dplyr::filter(term != "(Intercept)") |> 
  dplyr::select(x, term, estimate, std.error, p.value) |> 
  tidyr::pivot_wider(id_cols = x,
                     names_from = term,
                     values_from = c(estimate, std.error, p.value))

save(df, file = "out/predict_dbmi_baseline.Rdata")
load("out/predict_dbmi_baseline.Rdata")
df2 <- df |> dplyr::filter(!grepl("_clr", x)) |> tidyr::separate(x, "_", into = c('assay', 'var'), extra = 'merge')
df2 <- df2 |> dplyr::group_by(assay) |> 
  dplyr::mutate(padj = p.adjust(p.value_value, method = 'fdr'))



sig <- df2[df2$padj < 0.05 & !is.na(df2$padj),]
# Leptotrichia and defective NK Cells


# Plots:
load("data/data_raw.Rdata")
ASV_3jx_obf_72i
cells_single_cells_live_cd3neg_defective_freq_of_live
x <- as.data.frame(longFormat(data['cells_single_cells_live_cd3neg_defective_freq_of_live',,],
                              colDataCols = c('subjectId', 'interventionId', 'visitId', 'compliance', 'comprate', 'age_at_consent', 'bmi_at_consent', 'dbmi')))

x_m0 <- x |> dplyr::filter(visitId == 'M0' & interventionId != 'S')

library(ggplot2)

x_m0 |> dplyr::filter(!grepl("clr", assay)) |> 
  ggplot(aes(x = value,
             y = dbmi)) +
  geom_point() +
  ggpubr::stat_cor(method = 'spearman')


