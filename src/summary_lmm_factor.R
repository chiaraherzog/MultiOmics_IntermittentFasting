summary_lmm_factor <- function(dat, variables){
  
  library(parallel)
  
  cores <- parallel::detectCores()/4
  
  # Overall model ----------------
  cat('Overall model beginning ...')
  out <- mclapply(variables, function(i) {
    tryCatch(lmer(value ~ visitId*compliance + age_at_consent + (1 | subjectId),
                  data = dat[dat$rowname==i,],
                  REML = F),
             error = function(e) return(e))
  }, mc.cores = cores, mc.preschedule = T, mc.cleanup = T)
  names(out) <- variables
  
  ## Remove any with error
  out <- out |>
    purrr::keep( ~ !inherits(.x, 'error'))
  
  df <- dplyr::bind_rows(lapply(variables[variables %in% names(out)], function(x) as.data.frame(cbind(x, broom.mixed::tidy(out[[x]])))))
  
  df_all <- df |> 
    dplyr::filter(effect == 'fixed' & term != "(Intercept)") |> 
    dplyr::select(x, term, estimate, std.error, p.value) |> 
    tidyr::pivot_wider(id_cols = x,
                       names_from = term,
                       values_from = c(estimate, std.error, p.value))
  
  cat('done.\nHigh compliance only beginning ...')
  
  # High compliance only ----------------
  out <- mclapply(variables, function(i) {
    tryCatch(lmer(value ~ visitId + age_at_consent + (1 | subjectId),
                  data = dat[dat$rowname==i & dat$compliance=='high',],
                  REML = F),
             error = function(e) return(e))
  }, mc.cores = cores, mc.preschedule = T, mc.cleanup = T)
  names(out) <- variables
  
  ## Remove any with error
  out <- out |>
    purrr::keep( ~ !inherits(.x, 'error'))
  
  df_high <- dplyr::bind_rows(lapply(variables[variables %in% names(out)], function(x) as.data.frame(cbind(x, broom.mixed::tidy(out[[x]]))))) |> 
    dplyr::filter(effect == 'fixed' & term != "(Intercept)") |> 
    dplyr::select(x, term, estimate, std.error, p.value) |> 
    tidyr::pivot_wider(id_cols = x,
                       names_from = term,
                       values_from = c(estimate, std.error, p.value))
  
  # cat('done.\nLower compliance only beginning ...')
  
  # # Lower compliance only (for comparison) ----------------
  # out <- mclapply(variables, function(i) {
  #   tryCatch(lmer(value ~ time + age_at_consent + (1 | subjectId),
  #                 data = dat[dat$rowname==i & dat$compliance!='high',],
  #                 REML = F),
  #            error = function(e) return(e))
  # }, mc.cores = cores, mc.preschedule = T, mc.cleanup = T)
  # names(out) <- variables
  # 
  # ## Remove any with error
  # out <- out |>
  #   purrr::keep( ~ !inherits(.x, 'error'))
  # 
  # df_nhigh <- dplyr::bind_rows(lapply(variables[variables %in% names(out)], function(x) as.data.frame(cbind(x, broom.mixed::tidy(out[[x]]))))) |> 
  #   dplyr::filter(effect == 'fixed' & term != "Intercept") |> 
  #   dplyr::select(x, term, estimate, std.error, p.value) |> 
  #   tidyr::pivot_wider(id_cols = x,
  #                      names_from = term,
  #                      values_from = c(estimate, std.error, p.value))
  # 
  cat('done.\nIntervention interaction beginning ...')
  
  # Intervention ----------------
  out <- mclapply(variables, function(i) {
    tryCatch(lmer(value ~ visitId*interventionId + age_at_consent + (1 | subjectId),
                  data = dat[dat$rowname==i & dat$compliance == 'high',],
                  REML = F),
             error = function(e) return(e))
  }, mc.cores = cores, mc.preschedule = T, mc.cleanup = T)
  names(out) <- variables
  
  ## Remove any with error
  out <- out |>
    purrr::keep( ~ !inherits(.x, 'error'))
  
  df_intervent <- dplyr::bind_rows(lapply(variables[variables %in% names(out)], function(x) as.data.frame(cbind(x, broom.mixed::tidy(out[[x]]))))) |> 
    dplyr::filter(effect == 'fixed' & term != "(Intercept)") |> 
    dplyr::select(x, term, estimate, std.error, p.value) |> 
    tidyr::pivot_wider(id_cols = x,
                       names_from = term,
                       values_from = c(estimate, std.error, p.value))
  
  cat('done.\nIntervention sensitivity test beginning ...')
  
  # Intervention sensitivity test ----------------
  out <- mclapply(variables, function(i) {
    tryCatch(lmer(value ~ visitId*interventionId + age_at_consent + (1 | subjectId),
                  data = dat[dat$rowname==i & dat$compliance == 'high' & ((dat$interventionId == 'I') | (dat$interventionId == 'K' & dat$supplement == 'yes')),],
                  REML = F),
             error = function(e) return(e))
  }, mc.cores = cores, mc.preschedule = T, mc.cleanup = T)
  names(out) <- variables
  
  ## Remove any with error
  out <- out |>
    purrr::keep( ~ !inherits(.x, 'error'))
  
  df_intervent_sens <- dplyr::bind_rows(lapply(variables[variables %in% names(out)], function(x) as.data.frame(cbind(x, broom.mixed::tidy(out[[x]]))))) |> 
    dplyr::filter(effect == 'fixed' & term != "(Intercept)") |> 
    dplyr::select(x, term, estimate, std.error, p.value) |> 
    tidyr::pivot_wider(id_cols = x,
                       names_from = term,
                       values_from = c(estimate, std.error, p.value))
  
  cat('done.\nBMI included in basic model beginning ...')
  
  # Include BMI in the basic model ----------------
  out <- mclapply(variables, function(i) {
    tryCatch(lmer(value ~ visitId*compliance + age_at_consent + bmi_at_consent + (1 | subjectId),
                  data = dat[dat$rowname==i,],
                  REML = F),
             error = function(e) return(e))
  }, mc.cores = cores, mc.preschedule = T, mc.cleanup = T)
  names(out) <- variables
  
  ## Remove any with error
  out <- out |>
    purrr::keep( ~ !inherits(.x, 'error'))
  
  df <- dplyr::bind_rows(lapply(variables[variables %in% names(out)], function(x) as.data.frame(cbind(x, broom.mixed::tidy(out[[x]])))))
  
  df_bmi <- df |> 
    dplyr::filter(effect == 'fixed' & term != "(Intercept)") |> 
    dplyr::select(x, term, estimate, std.error, p.value) |> 
    tidyr::pivot_wider(id_cols = x,
                       names_from = term,
                       values_from = c(estimate, std.error, p.value))
  
  # Full model --------------
  cat('Full model starting...')

  out <- mclapply(variables, function(i) {
    tryCatch(lmer(value ~ age_at_consent + bmi_at_consent + visitId*compliance*interventionId + (1|subjectId),
                  data = dat[dat$rowname==i,],
                  REML = F),
             error = function(e) return(e))
  }, mc.cores = cores, mc.preschedule = T, mc.cleanup = T)
  names(out) <- variables
  
  ## Remove any with error
  out <- out |>
    purrr::keep( ~ !inherits(.x, 'error'))
  
  df <- dplyr::bind_rows(lapply(variables[variables %in% names(out)], function(x) as.data.frame(cbind(x, broom.mixed::tidy(out[[x]])))))
  
  df_full <- df |> 
    dplyr::filter(effect == 'fixed' & term != "(Intercept)") |> 
    dplyr::select(x, term, estimate, std.error, p.value) |> 
    tidyr::pivot_wider(id_cols = x,
                       names_from = term,
                       values_from = c(estimate, std.error, p.value))
  
  cat('done.\nExporting.')
  
  out_lmm_factor <- list(overall_interaction = df_all,
                  high_compliance = df_high,
                  # lower_compliance = df_nhigh,
                  intervention = df_intervent,
                  intervent_sens = df_intervent_sens,
                  bmi_sens = df_bmi,
                  full = df_full)
  
  save(out_lmm_factor, file = "out/lmm_factor.Rdata")
  
  
  return(out_lmm)
}
