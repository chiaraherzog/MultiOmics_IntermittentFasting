load("data/data_raw.Rdata")
df_raw <- as.data.frame(longFormat(data[vars,,],
                                   colData = c("age_at_consent", "interventionId", "subjectId", "visitId", "time", 'compliance', 'comprate', 'mpstatrs', 'pattern', 'dropout_date'))) |> 
  dplyr::filter(!is.na(value))

df_raw <- df_raw |>
  dplyr::filter(interventionId!='S' & !visitId %in% c("M12", "M18")) |>
  dplyr::mutate(compliance = factor(compliance, levels = c("low", 'medium', 'high')))

# Load data (change, raw)
load("data/data_baseline_change.Rdata")
df_change <- as.data.frame(longFormat(data[vars,,], colData = c("age_at_consent", "interventionId", "subjectId", "visitId", "time", 'compliance', 'comprate', 'pattern', 'mpstatrs', 'dropout_date')))

df_change <- df_change |>
  dplyr::filter(interventionId!='S' & !visitId %in% c("M12", "M18")) |>
  dplyr::mutate(compliance = factor(compliance, levels = c("low", 'medium', 'high')))|> 
  dplyr::filter(!is.na(value))


df_raw <- df_raw |> dplyr::mutate(dropout = ifelse(is.na(dropout_date), "no", "yes"))
df_change <- df_change |> dplyr::mutate(dropout = ifelse(is.na(dropout_date), "no", "yes"))

df_change |>dplyr::filter(rowname == 'fm') |> 
  ggplot(aes(x = dropout,
                               y = value))+
  geom_boxplot() +
  facet_wrap(~visitId) +
  ggpubr::stat_compare_means(ref.group = 'no')

df_raw |>dplyr::filter(rowname == 'weight') |> q
  ggplot(aes(x = dropout,
             y = value))+
  geom_boxplot() +
  facet_wrap(~visitId) +
  ggpubr::stat_compare_means(ref.group = 'no')
