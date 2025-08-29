library(dplyr)
library(data.table)
## wearable
load("~/Dropbox/data/tirolgesund/wearable/wearable_activity.Rdata")
wearable <- data

## mae
library(MultiAssayExperiment)
load("data/data_raw.Rdata")
pheno <- colData(data)
t <- data@metadata
dates <- as.data.frame(t$`timing of collection`) |> dplyr::filter(assay == 'Blood haemogram') |> 
  dplyr::mutate(subjectId = substr(primary, 1, 4),
                visitId = substr(primary, 5, nchar(primary))) |> 
  dplyr::rename(visit_day_t = t)

# make sure your timestamps are Dates (or POSIXct that can be compared as dates)
# assuming wearable$t and dates$visit_day_t are Date; if not, cast with as.Date()

# 1) Build 14-day windows ---------------------------------

# Baseline (M0): first 14 days of wearable per subject
baseline_windows <- wearable |> 
  group_by(subjectId) |> 
  summarise(
    start = 1,
    end   = 14,
    label = "M0",
    .groups = "drop"
  )

# Visit windows: 14 days before each visit (exclude M0 so baseline is “first 2 weeks”)
visit_windows <- dates |> 
  filter(!is.na(visit_day_t), !grepl("^M0$", visitId, ignore.case = TRUE)) |> 
  transmute(
    subjectId,
    start = visit_day_t - 13,
    end   = visit_day_t,
    label = visitId
  )

interval_windows <- bind_rows(baseline_windows, visit_windows)

# 2) Assign each daily wearable row to a window ------------
nonwear_threshold <- 500
active_threshold  <- 7000

wearable_filtered <- wearable |> 
  filter(steps > nonwear_threshold) |> 
  mutate(activity_level = ifelse(steps >= active_threshold, "Active", "Sedentary"))

d_int  <- as.data.table(interval_windows)
d_days <- as.data.table(wearable_filtered)

setkey(d_int, subjectId, start, end)

# t in (start, end] — i.e., strictly after start, up to and including end.
# (Prevents double-assign at boundaries if two windows abut.)
d_days[d_int, on = .(subjectId, t > start, t <= end), label := i.label]

dailies <- d_days[!is.na(label)]


interval_summary <- dailies |> 
  group_by(subjectId, label) |> 
  summarise(
    across(c(activekcal, bmr, fitnessage, remsleep, rhr, steps, minvigint),
           ~ mean(.x, na.rm = TRUE),
           .names = "mean_{.col}"),
    activity_level_summary = {
      x <- activity_level[!is.na(activity_level)]
      if (length(x) == 0) NA_character_
      else {
        tab <- table(x)
        # pick the mode; if tied, keep the first seen in the original order
        modes <- names(tab)[tab == max(tab)]
        x[match(TRUE, x %in% modes)]
      }
    },
    .groups = "drop"
  )

# append pheno
interval_summary <- interval_summary |>
  dplyr::rename(visitId = label) |> dplyr::left_join(pheno, copy = T)

# Individuals who remained in the study (no dropout, complete data for steps)
complete <- interval_summary |> 
  dplyr::filter(is.na(dropout_reason)) |> 
  dplyr::group_by(subjectId) |> 
  dplyr::filter(all(c("M0", "M2", "M4", "M6") %in% visitId)) |>
  dplyr::pull(subjectId) |> 
  unique()

interval_summary |> 
  dplyr::filter(interventionId != 'S') |> 
  dplyr::filter(subjectId %in% complete) |> 
  ggplot(aes(x = visitId,
             y = mean_steps)) +
  geom_boxplot(outlier.shape = NA) +
  geom_line(aes(group = subjectId),
            alpha = 0.6,
            colour = 'grey50') +
  facet_wrap(~interventionId) +
  ggpubr::stat_compare_means(ref.group = 'M0', paired = T,
                             label = 'p.format',label.y.npc = 0.95) +
  theme_bw() +
  labs(y = 'average steps/day\nduring interval')



interval_summary |> 
  dplyr::filter(interventionId != 'S') |> 
  dplyr::mutate(compliance = factor(compliance, levels = c("low", "medium", "high"))) |> 
  dplyr::filter(subjectId %in% complete) |> 
  ggplot(aes(x = visitId,
             y = mean_steps,
             fill = compliance)) +
  geom_boxplot(outlier.shape = NA,
               alpha = 0.4) +
  geom_point(aes(colour = compliance),
             position = position_dodge(width = 0.75),
             size = 1,
             alpha = 0.8) +
  # geom_line(aes(group = subjectId),
  #           alpha = 0.6,
  #           colour = 'grey50') +
  facet_wrap(~interventionId) +
  theme_bw() +
  labs(y = 'average steps/day\nduring interval')


act <- interval_summary |> 
  dplyr::filter(interventionId != 'S') |> 
  dplyr::filter(is.na(dropout_reason)) |> 
  dplyr::group_by(visitId) |> 
  dplyr::mutate(n_total = n()) |> 
  dplyr::group_by(visitId, interventionId, activity_level_summary) |> 
  dplyr::reframe(n = n(),
                 prop = n/n_total) |> dplyr::distinct()


act |> 
  ggplot(aes(x = visitId,
             y = prop,
             fill = activity_level_summary)) +
  geom_col()

ggplot(aes(x = visitId,
           y = mean_minvigint,
           fill  = compliance_smkgroup)) +
  geom_boxplot() +
  facet_wrap(~interventionId)


complete <- interval_summary |> 
  dplyr::filter(!is.na(mean_rhr)) |> 
  dplyr::filter(is.na(dropout_reason)) |> 
  dplyr::group_by(subjectId) |> 
  dplyr::filter(all(c("M0", "M2", "M4", "M6") %in% visitId)) |>
  dplyr::pull(subjectId) |> 
  unique()

interval_summary |> 
  dplyr::filter(interventionId != 'S') |> 
  dplyr::filter(subjectId %in% complete) |> 
  ggplot(aes(x = visitId,
             y = mean_rhr)) +
  geom_boxplot(outlier.shape = NA) +
  geom_line(aes(group = subjectId),
            alpha = 0.6,
            colour = 'grey50') +
  facet_wrap(~interventionId) +
  ggpubr::stat_compare_means(ref.group = 'M0', paired = T,
                             label = 'p.format',label.y.npc = 0.95) +
  theme_bw() +
  labs(y = 'average steps/day\nduring interval')

interval_summary |> 
  dplyr::filter(interventionId != 'S') |> 
  dplyr::mutate(compliance = factor(compliance, levels = c("low", "medium", "high"))) |> 
  dplyr::filter(subjectId %in% complete) |> 
  ggplot(aes(x = visitId,
             y = mean_rhr,
             fill = compliance)) +
  geom_boxplot(outlier.shape = NA,
               alpha = 0.4) +
  geom_point(aes(colour = compliance),
             position = position_dodge(width = 0.75),
             size = 1,
             alpha = 0.8) +
  # geom_line(aes(group = subjectId),
  #           alpha = 0.6,
  #           colour = 'grey50') +
  # facet_wrap(~compliance) +
  theme_bw() +
  labs(y = 'average steps/day\nduring interval')


complete <- interval_summary |> 
  dplyr::filter(!is.na(mean_remsleep)) |> 
  dplyr::filter(is.na(dropout_reason)) |> 
  dplyr::group_by(subjectId) |> 
  dplyr::filter(all(c("M0", "M2", "M4", "M6") %in% visitId)) |>
  dplyr::pull(subjectId) |> 
  unique()

interval_summary |> 
  dplyr::filter(interventionId != 'S') |> 
  dplyr::filter(subjectId %in% complete) |> 
  ggplot(aes(x = visitId,
             y = mean_minvigint)) +
  geom_boxplot(outlier.shape = NA) +
  geom_line(aes(group = subjectId),
            alpha = 0.6,
            colour = 'grey50') +
  # facet_wrap(~interventionId) +
  ggpubr::stat_compare_means(ref.group = 'M0', paired = T,
                             label = 'p.format',label.y.npc = 0.95) +
  theme_bw() +
  labs(y = 'average steps/day\nduring interval')

interval_summary |> 
  dplyr::filter(interventionId != 'S') |> 
  dplyr::mutate(compliance = factor(compliance, levels = c("low", "medium", "high"))) |> 
  dplyr::filter(subjectId %in% complete) |> 
  ggplot(aes(x = visitId,
             y = mean_minvigint,
             fill = compliance)) +
  geom_boxplot(outlier.shape = NA,
               alpha = 0.4) +
  geom_point(aes(colour = compliance),
             position = position_dodge(width = 0.75),
             size = 1,
             alpha = 0.8) +
  # geom_line(aes(group = subjectId),
  #           alpha = 0.6,
  #           colour = 'grey50') +
  # facet_wrap(~compliance) +
  theme_bw() +
  labs(y = 'average steps/day\nduring interval')
