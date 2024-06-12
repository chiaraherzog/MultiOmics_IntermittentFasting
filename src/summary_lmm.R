summary_lmm <- function(dat, variables){

  library(parallel)
  
  cores <- floor(parallel::detectCores()/4)
  
  # Overall model ----------------
  cat('Overall model beginning ...')
  out <- mclapply(variables, function(i) {
    tryCatch(lmer(value ~ time*compliance + age_at_consent + (1 | subjectId),
                  data = dat[dat$rowname==i,],
                  REML = F),
             error = function(e) return(e))
  }, mc.cores = cores, mc.preschedule = T, mc.cleanup = T)
  names(out) <- variables
  
  ## Remove any with error
  out <- out |>
    purrr::keep( ~ !inherits(.x, 'error'))
  
  v <- variables[variables %in% names(out)]
  
  # Bind
  df <- dplyr::bind_rows(lapply(v, function(x) as.data.frame(cbind(x, broom.mixed::tidy(out[[x]])))))
  
  df_all <- df |> 
    dplyr::filter(effect == 'fixed' & term %in% c("time", "age_at_consent", "time:compliancehigh")) |> 
    dplyr::select(x, term, estimate, std.error, p.value) |> 
    tidyr::pivot_wider(id_cols = x,
                       names_from = term,
                       values_from = c(estimate, std.error, p.value))
  
  # High compliance only ----------------
  cat('done.\nHigh compliance only beginning ...')
  
  out <- mclapply(variables, function(i) {
    tryCatch(lmer(value ~ time + age_at_consent + (1 | subjectId),
                  data = dat[dat$rowname==i & dat$compliance=='high',],
                  REML = F),
             error = function(e) return(e))
  }, mc.cores = cores, mc.preschedule = T, mc.cleanup = T)
  names(out) <- variables
  
  ## Remove any with error
  out <- out |>
    purrr::keep( ~ !inherits(.x, 'error'))
  
  v <- variables[variables %in% names(out)]
    
  # Bind
  df_high<- dplyr::bind_rows(lapply(v, function(x) as.data.frame(cbind(x, broom.mixed::tidy(out[[x]]))))) |> 
    dplyr::filter(effect == 'fixed' & term %in% c("time", "age_at_consent")) |> 
    dplyr::select(x, term, estimate, std.error, p.value) |> 
    tidyr::pivot_wider(id_cols = x,
                       names_from = term,
                       values_from = c(estimate, std.error, p.value))
  
  
  # Lower compliance only (for comparison) ----------------
  cat('done.\nLower compliance only beginning ...')
  out <- mclapply(variables, function(i) {
    tryCatch(lmer(value ~ time + age_at_consent + (1 | subjectId),
                  data = dat[dat$rowname==i & dat$compliance!='high',],
                  REML = F),
             error = function(e) return(e))
  }, mc.cores = cores, mc.preschedule = T, mc.cleanup = T)
  names(out) <- variables
  
  ## Remove any with error
  out <- out |>
    purrr::keep( ~ !inherits(.x, 'error'))
  
  v <- variables[variables %in% names(out)]
  
  ## Bind
  df_nhigh <- dplyr::bind_rows(lapply(v, function(x) as.data.frame(cbind(x, broom.mixed::tidy(out[[x]]))))) |> 
    dplyr::filter(effect == 'fixed' & term %in% c("time", "age_at_consent")) |> 
    dplyr::select(x, term, estimate, std.error, p.value) |> 
    tidyr::pivot_wider(id_cols = x,
                       names_from = term,
                       values_from = c(estimate, std.error, p.value))

  
  # Intervention ---------------- 
  cat('done.\nIntervention interaction beginning ...')
  
  out <- mclapply(variables, function(i) {
    tryCatch(lmer(value ~ time*interventionId + age_at_consent + (1 | subjectId),
                  data = dat[dat$rowname==i & dat$compliance == 'high',],
                  REML = F),
             error = function(e) return(e))
  }, mc.cores = cores, mc.preschedule = T, mc.cleanup = T)
  names(out) <- variables
  
  ## Remove any with error
  out <- out |>
    purrr::keep( ~ !inherits(.x, 'error'))
  
  v <- variables[variables %in% names(out)]
  
  df_intervent <- dplyr::bind_rows(lapply(v, function(x) as.data.frame(cbind(x, broom.mixed::tidy(out[[x]]))))) |> 
    dplyr::filter(effect == 'fixed' & term %in% c("time", "age_at_consent", "time:interventionIdK")) |> 
    dplyr::select(x, term, estimate, std.error, p.value) |> 
    tidyr::pivot_wider(id_cols = x,
                       names_from = term,
                       values_from = c(estimate, std.error, p.value))
  
  # Intervention supplement sensitivity test ---------------- 
  cat('done.\nIntervention sensitivity test beginning ...')
  out <- mclapply(variables, function(i) {
    tryCatch(lmer(value ~ time*interventionId + age_at_consent + (1 | subjectId),
                  data = dat[dat$rowname==i & dat$compliance == 'high' & ((dat$interventionId == 'I') | (dat$interventionId == 'K' & dat$supplement == 'yes')),],
                  REML = F),
             error = function(e) return(e))
  }, mc.cores = cores, mc.preschedule = T, mc.cleanup = T)
  names(out) <- variables
  
  ## Remove any with error
  out <- out |>
    purrr::keep( ~ !inherits(.x, 'error'))
  
  v <- variables[variables %in% names(out)]
  
  df_intervent_sens <- dplyr::bind_rows(lapply(v, function(x) as.data.frame(cbind(x, broom.mixed::tidy(out[[x]]))))) |> 
    dplyr::filter(effect == 'fixed' & term %in% c("time", "age_at_consent", "time:interventionIdK")) |> 
    dplyr::select(x, term, estimate, std.error, p.value) |> 
    tidyr::pivot_wider(id_cols = x,
                       names_from = term,
                       values_from = c(estimate, std.error, p.value))
  
  cat('done.\nBMI included in basic model beginning ...')
  
  # Overall model with BMI ----------------
  out <- mclapply(variables, function(i) {
    tryCatch(lmer(value ~ time*compliance + age_at_consent + bmi_at_consent + (1 | subjectId),
                  data = dat[dat$rowname==i,],
                  REML = F),
             error = function(e) return(e))
  }, mc.cores = cores, mc.preschedule = T, mc.cleanup = T)
  names(out) <- variables
  
  ## Remove any with error
  out <- out |>
    purrr::keep( ~ !inherits(.x, 'error'))
  
  v <- variables[variables %in% names(out)]
  
  # Bind
  df <- dplyr::bind_rows(lapply(v, function(x) as.data.frame(cbind(x, broom.mixed::tidy(out[[x]])))))
  
  df_bmi <- df |> 
    dplyr::filter(effect == 'fixed' & term %in% c("time", "age_at_consent", "time:compliancehigh")) |> 
    dplyr::select(x, term, estimate, std.error, p.value) |> 
    tidyr::pivot_wider(id_cols = x,
                       names_from = term,
                       values_from = c(estimate, std.error, p.value))
  
  
  cat('done.\nExporting.')
  
  out_lmm <- list(overall_interaction = df_all,
                  high_compliance = df_high,
                  lower_compliance = df_nhigh,
                  intervention = df_intervent,
                  intervention_sensitivity = df_intervent_sens,
                  overall_bmi = df_bmi)
  return(out_lmm)
}
