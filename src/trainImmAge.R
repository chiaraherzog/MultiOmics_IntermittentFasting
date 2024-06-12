trainImmAge <- function(data, exp, populations,
                        mainOnly = F){
  
  # Extract data ------
  tr <- as.data.frame(wideFormat(data[,data@colData$visitId == 'M0' & data@colData$interventionId!='S',
                                      exp],
                                 colData = c('age_at_consent'))) |> tidyr::drop_na()
  colnames(tr) <- gsub("Flow.cytometry..T.cell.staining_|Flow.cytometry..white.blood.cell.staining_", "", colnames(tr))
  
  trdata <- tr |> dplyr::select(-c('age_at_consent', 'primary'))
  
  if(mainOnly == T){
    trdata <- trdata |> 
      dplyr::select(any_of(dplyr::filter(populations, main_analysis == 'yes')$name))
  }
  
  trpheno <- tr |> dplyr::select(c('age_at_consent'))
  
  # We train the immAge using the agePred function and select maximum correlation with age: ------
  models <- list("ridge" = 0,
                 "elnet" = 0.5,
                 "lasso" = 1)
  
  res <- lapply(models, function(x){
    agePred(data_tr = t(trdata),
            age_tr = as.numeric(trpheno$age_at_consent),
            data_val = t(trdata),
            age_val = as.numeric(trpheno$age_at_consent),
            alpha = x)
  })
  
  cors <- lapply(names(res), function(m){
    as.numeric(cor(trpheno$age_at_consent, res[[m]]$val_predictor))
  })

  keep <- which.max(cors)
  cat(names(res)[keep], ' model has the largest correlation.\n', sep = '')
  
  res <- res[[keep]]
  intercept <- as.numeric(coef(res$fit.cv, s = "lambda.min")[1]) 
  
  coefs <- as.matrix(coef(res$fit.cv, s = "lambda.min")) |> 
    as.data.frame() |> 
    tibble::rownames_to_column('name') |> 
    dplyr::filter(s1 != 0 & name != '(Intercept)') |> 
    dplyr::arrange(desc(abs(s1)))
  
  dt <- DT::datatable(coefs |> 
                  dplyr::left_join(dplyr::select(populations, name, `population name`)) |> 
                  dplyr::select(`population name`, s1))

  # Compute values in the full dataset -------
  full <- as.data.frame(wideFormat(data[,,
                                        exp],
                                        colData = c('interventionId', 'subjectId', 'visitId', 'compliance', 'age_at_consent', 'bmi_at_consent', 'smoking_py')))
  
  colnames(full) <- gsub("Flow.cytometry..T.cell.staining_|Flow.cytometry..white.blood.cell.staining_", "", colnames(full))
  
  full <- full |>
    dplyr::mutate(age = case_match(visitId,
                                   "M0" ~ age_at_consent,
                                   "M2" ~ age_at_consent + 2/12,
                                   "M4" ~ age_at_consent + 4/12,
                                   "M6" ~ age_at_consent + 0.5))
  
  
  dat_fu <- full |> 
    dplyr::select(coefs$name) |> 
    as.matrix()
  
  if(!identical(colnames(dat_fu), coefs$name)){
    stop("colnames of full dataset not identical with coefs. please investigate.\n")
  }
  
  full$ImmAge <- apply(dat_fu, 1, function(x){
    intercept + sum(x * coefs$s1)
  })
  
  full <- full |> 
    dplyr::rowwise() |> 
    dplyr::mutate(AgeDev = ImmAge - age) |> 
    dplyr::ungroup()
  
  # Initial plots: ---------
  corr_age <- full |> 
    ggplot(aes(x = age,
               y = ImmAge)) +
    geom_point() +
    geom_smooth(method = 'lm') +
    facet_wrap(~visitId,
               nrow = 1) +
    ggpubr::stat_cor() +
    labs(subtitle = 'Correlation with age')
  
  corr_age_intervent <- full |> 
    dplyr::mutate(study_arm = ifelse(interventionId == 'S', 'S', 'IF')) |> 
    ggplot(aes(x = age,
               y = ImmAge,
               colour = study_arm)) +
    geom_point() +
    geom_smooth(method = 'lm') +
    facet_wrap(~visitId,
               nrow = 1) +
    ggpubr::stat_cor() +
    labs(subtitle = 'Correlation with age by intervention')
  
  corr_agedev_intervent <- full |> 
    dplyr::mutate(study_arm = ifelse(interventionId == 'S', 'S', 'IF')) |> 
    ggplot(aes(x = age,
               y = AgeDev,
               colour = study_arm)) +
    geom_point() +
    geom_smooth(method = 'lm') +
    facet_wrap(~visitId,
               nrow = 1) +
    ggpubr::stat_cor() +
    labs(subtitle = 'Correlation with AgeDev by intervention')
  
  fit <- lm(ImmAge ~ age_at_consent, data = full[full$visitId=='M0',])
  full$ImmAge_adj <- full$ImmAge - predict(fit, full)
  
  fit <- lm(AgeDev ~ age_at_consent, data = full[full$visitId=='M0',])
  full$AgeDev_adj <- full$AgeDev - predict(fit, full)
  
  bmi <- full |> 
    dplyr::filter(visitId == 'M0') |> 
    ggplot(aes(x = bmi_at_consent,
               y = ImmAge_adj)) +
    geom_point() +
    geom_smooth(method = 'lm') +
    facet_wrap(~visitId,
               nrow = 1) +
    ggpubr::stat_cor() +
    labs(subtitle = 'Correlation of AgeDev with BMI (M0)')
  
  spy <- full |> 
    dplyr::filter(visitId == 'M0') |> 
    ggplot(aes(x = smoking_py,
               y = ImmAge_adj)) +
    geom_point() +
    geom_smooth(method = 'lm') +
    facet_wrap(~visitId,
               nrow = 1) +
    ggpubr::stat_cor() +
    labs(subtitle = 'Correlation with cigarettes (M0)')
  

  adj <- full |> 
    dplyr::filter(visitId == 'M0') |> 
    dplyr::mutate(study_arm = ifelse(interventionId == 'S', 'S', 'IF')) |> 
    ggplot(aes(x = study_arm,
               y = ImmAge_adj)) +
    geom_boxplot() +
    ggpubr::stat_compare_means(comparisons = list(c("IF", "S")),
                               label = 'p.format') +
    labs(subtitle = 'M0 values')
  
# Follow-up visits (adjusted)
  follow_raw <- full |> 
    dplyr::mutate(compliance = factor(compliance, levels = c("low", 'medium', 'high', 'lower compliance', 'higher compliance'))) |> 
    ggplot(aes(x = visitId,
               y = ImmAge_adj)) +
    geom_boxplot() +
    geom_line(aes(group = subjectId),
              alpha = 0.1) +
    facet_wrap(interventionId~compliance) +
    ggpubr::stat_compare_means(ref.group = 'M0',
                               label = 'p.signif',
                               label.y.npc = 0.95) +
    labs(subtitle = 'Raw values over time')
  
  
  follow_paired <- full |> 
    dplyr::mutate(compliance = factor(compliance, levels = c("low", 'medium', 'high', 'lower compliance', 'higher compliance'))) |> 
    tidyr::pivot_wider(id_cols = c(subjectId, compliance, age_at_consent, interventionId),
                       names_from = visitId,
                       values_from = ImmAge_adj) |> 
    dplyr::rowwise() |> 
    dplyr::mutate(M2 = M2-M0,
                  M4 = M4-M0,
                  M6 = M6-M0,
                  M0 = 0) |>
    dplyr::ungroup() |> 
    dplyr::filter(!is.na(M0) & !is.na(M2) & !is.na(M6) & !is.na(M4)) |> 
    tidyr::pivot_longer(cols = M0:M6,
                        names_to = 'visitId',
                        values_to = 'ImmAge') |> 
    ggplot(aes(x = visitId,
               y = ImmAge)) +
    geom_boxplot() +
    geom_line(aes(group = subjectId),
              alpha = 0.1) +
    facet_wrap(interventionId ~ compliance) +
    ggpubr::stat_compare_means(ref.group = 'M0',
                               label = 'p.signif', paired = T,
                               label.y.npc = 0.95) +
    labs(y = '∆ ImmAge (adjusted)\nfrom baseline')
  
  
  out <- list(data = full,
              res = res,
              coef = dt,
              plots = list("corr_age" = corr_age,
                           "corr_age_intervent" = corr_age_intervent,
                           "corr_agedev_intervent" = corr_agedev_intervent,
                           "bmi" = bmi,
                           "spy" = spy,
                           "adj" = adj,
                           "follow_raw" = follow_raw,
                           "follow_paired" = follow_paired))
  
  return(out)
}
