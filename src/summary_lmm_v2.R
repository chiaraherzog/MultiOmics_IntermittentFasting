#' @param dat dataframe to input
#' @param variables variables to loop over
#' @param models list of models (optional). If no list supplied as per default, standard models are run.
#' @param timeType should time be treated as continuous or as a factor?
#' @param outPath path to place the output of files
#' @param outName output fie name. default is 'lmm.Rdata'
#' @param cores cores to allocate for parallel processing. NULL by default, which means 1/4 of all available cores are used by default.
#' @returns list of LMM results for each model, within the list models are supplied as dataframes.

summary_lmm <- function(dat,
                        variables,
                        models = NULL,
                        timeType = NULL,
                        outPath = file.path(getwd(), "out"),
                        outName = "lmm.Rdata",
                        cores = NULL
                        ){
  
  # Libraries
  library(parallel)
  library(dplyr)
  library(broom.mixed)
  library(R.utils)
  
  # Set up cores
  if(is.null(cores)){
  cores <- parallel::detectCores()/8
  }
  
  cat("Starting. You are using **", cores, "** cores and are looping over ", length(variables), " variables.\n", sep = "")
  
  # Set up variables
  if(is.null(models) & is.null(timeType)){
    stop("[SETUP ERROR!]\nNo models or default time type defined.\nPlease provide either a list of models to run or set timeType to continuous or factor.")
  }
  
  if(is.null(models) & timeType == 'continuous'){
    models <- list("Minimal model" = 'value ~ time + age_at_consent + bmi_at_consent + (1|subjectId)',
                   "Basic model" = 'value ~ time*compliance + age_at_consent + (1 | subjectId)',
                   "Basic model with BMI" = 'value ~ time*compliance + age_at_consent + bmi_at_consent + (1 | subjectId)',
                   "Basic model with BMI and BMI change" = 'value ~ time*compliance + age_at_consent + dbmi + (1 | subjectId)',
                   "Intervention (higher compliance)" = 'value ~ time*interventionId + age_at_consent + (1 | subjectId)',
                   "Intervention (higher compliance, supplement sensitivity test)" = 'value ~ time*interventionId + age_at_consent + (1 | subjectId)',
                   "Full model" = 'value ~ age_at_consent + bmi_at_consent + time*compliance*interventionId + (1|subjectId)'
    )
  } else if (is.null(models) & timeType == 'factor'){
    models <- list("Minimal model" = 'value ~ visitId + age_at_consent + bmi_at_consent + (1|subjectId)',
                   "Minimal model (ordered)" = 'value ~ visitIdOrdered + age_at_consent + bmi_at_consent + (1|subjectId)',
                   
                   # Basic models
                   "Basic model (no BMI)" = 'value ~ visitId*compliance + age_at_consent + (1 | subjectId)',
                   "Basic model with BMI" = 'value ~ visitId*compliance + age_at_consent + bmi_at_consent + (1 | subjectId)',
                   "Basic model with BMI and BMI change" = 'value ~ visitId*compliance + age_at_consent + bmi_at_consent + dbmi + (1 | subjectId)',
                   "Basic model with BMI - ordered" = 'value ~ visitIdOrdered*compliance + age_at_consent + bmi_at_consent + (1 | subjectId)',
                   "Basic model with BMI and BMI change - ordered" = 'value ~ visitIdOrdered*compliance + age_at_consent + bmi_at_consent + dbmi + (1 | subjectId)',
                   
                   # Intervention
                   "Intervention (higher compliance)" = 'value ~ visitId*interventionId + age_at_consent + bmi_at_consent + (1 | subjectId)',
                   "Intervention (higher compliance supplement sensitivity test)" = 'value ~ visitId*interventionId + age_at_consent + bmi_at_consent + (1 | subjectId)',
                   "Intervention (all)" = 'value ~ visitId*interventionId + age_at_consent + bmi_at_consent + (1|subjectId)',
                   "Intervention (higher compliance) - ordered" = 'value ~ visitIdOrdered*interventionId + age_at_consent + bmi_at_consent + (1 | subjectId)',
                   
                   "Full model" = 'value ~ age_at_consent + bmi_at_consent + visitId*compliance*interventionId + (1|subjectId)',
                   
                   "Menopause interaction" = 'value ~ age_at_consent + bmi_at_consent + visitId*mpstatrs + interventionId + (1|subjectId)',
                   "Menopause general" = "value ~ age_at_consent + bmi_at_consent + visitId + mpstatrs + interventionId + (1|subjectId)",
                   "Basic model with BMI, Menopause and Age" = 'value ~ visitId*compliance + age_at_consent + bmi_at_consent + mpstatrs + (1 | subjectId)'
    )
  } else if(!is.null(models)){
    models <- models
  }
  
  cat("Models defined.\n\n")
  
  sink(file.path(outPath, paste0(sub("\\..*$", "", outName), "_models.txt")))
  print(models)
  sink()
  
  
  # Set up outList and loop over models
  
  outList <- vector("list", length(models))
  
  for (m in 1:length(models)){
    cat(names(models)[m], "beginning ... ")
    
    # Subset right population depending on model requirements; needs more work to accommodate custom models.
    if(grepl("Basic|Full|Minimal|all", names(models)[m])){
      tmp <- dat
    } else if (grepl("supplement", names(models)[m])){
      tmp <- dat |> 
        dplyr::filter(compliance == 'high' & ((interventionId == 'I') | (interventionId == 'K' & supplement == 'yes')))
    } else if (grepl("higher", names(models)[m]) | grepl("Menopause", names(models)[m])){
      tmp <- dat |> 
        dplyr::filter(compliance == 'high')
    }
    
    out <- mclapply(variables, function(i) {
      tryCatch(lmer(models[[m]],
                    data = tmp[tmp$rowname==i,],
                    REML = F),
               error = function(e) return(e))
    }, mc.cores = cores, mc.preschedule = T, mc.cleanup = T)
    names(out) <- variables
    
    ## Remove any with error
    out <- out |>
      purrr::keep( ~ !inherits(.x, 'error')) 
    
    ## Bind
    
    df <- dplyr::bind_rows(lapply(variables[variables %in% names(out)], function(x) as.data.frame(cbind(x, broom.mixed::tidy(out[[x]]))))) |> 
      dplyr::filter(effect == 'fixed' & term != "(Intercept)") |> 
      dplyr::select(x, term, estimate, std.error, p.value) |> 
      tidyr::pivot_wider(id_cols = x,
                         names_from = term,
                         values_from = c(estimate, std.error, p.value))
    
    ## Append to outList
    outList[[m]] <- df
    
    # Save interim dataframe
    save(df, file = file.path(outPath, paste0(R.utils::toCamelCase(names(models)[m]), "_temp.Rdata")))
    
    cat("done.\n")
    
  }
  
  cat('done.\nExporting.')
  
  out_lmm <- outList
  names(out_lmm) <- names(models)
  save(out_lmm, file = file.path(outPath, outName))
  
  # Delete intermediary files
  for (n in names(models)){
    file.remove(file.path(outPath, paste0(R.utils::toCamelCase(n), "_temp.Rdata")))
  }
  
  return(out_lmm)
}
