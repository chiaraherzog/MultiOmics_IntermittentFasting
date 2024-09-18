#' Permutation testing of Longitudinal Data
#'
#' This function performs a permutation test on a longitudinal dataset by shuffling the `visitId` within each `subjectId` and `rowname` group,
#' fitting a specified mixed-effects model for each unique `rowname`, and repeating this process a specified number of times.
#' The results from each iteration are combined into a single data frame.
#'
#' @param df A data frame containing the longitudinal data, with columns for `subjectId`, `rowname`, `visitId`, and other variables required by the model.
#' @param n An integer specifying the number of iterations to perform for the permutation test. Default is 500.
#' @param model A string representing the formula for the mixed-effects model. The default is `'value ~ visitId*compliance + age_at_consent + bmi_at_consent + (1 | subjectId)'`.
#' @param variables A character vector specifying the subset of `rowname` values to include in the analysis. If NULL (default), all `rowname` values are used.
#' @return A data frame containing the results of the permutation test, including estimates, standard errors, and p-values for each term in the model,
#'         for each `rowname`, across all iterations. The data frame also includes a `replicate_id` column to indicate the iteration number.
#' @examples
#' \dontrun{
#'   # Assuming you have a data frame `df` with the required columns:
#'   results <- cleaned_shuffle_and_model(df, n = 500, model = 'value ~ visitId*compliance + age_at_consent + bmi_at_consent + (1 | subjectId)')
#'   print(results)
#' }
#' @import dplyr data.table future progressr broom.mixed tidyr purrr lme4 cli
#' @export

permutationTesting <- function(df,
                               variables = NULL,
                               n = 500,
                               model = 'value ~ visitId*compliance + age_at_consent + bmi_at_consent + (1 | subjectId)') {
  
  library(dplyr)
  library(data.table)
  library(future)
  library(future.apply)
  library(progressr)
  library(broom.mixed)
  library(tidyr)
  library(purrr)
  library(lme4)
  library(cli)
  library(parallel)
  
  set.seed(1234)
  
  # Set up parallel processing plan
  plan(multisession, workers = availableCores()/8)
  
  # Set up progress bar handler
  handlers(handler_txtprogressbar(char = cli::col_red(cli::symbol$heart), width = 50))
  
  # Subset the data to include only the specified rownames if provided
  if (!is.null(variables)) {
    df <- df[df$rowname %in% variables, ]
  }
  
  shuffle_and_model_iteration <- function(df, model, variables) {
    
    # Shuffle 'visitId' within each 'subjectId' and 'rowname' group; leveraging data.table package for efficiency
    setDT(df)
    
    if(!grepl("interventionId", model)){
    df[, visitId := if (.N > 1) sample(visitId) else visitId, by = .(subjectId, rowname)]
    } else {
      df[, interventionId := if (.N > 1) {
        # Assign either "I" or "K" to all rows within the group
        chosen_intervention <- sample(c("I", "K"), 1)
        chosen_intervention
      } else {
        interventionId
      }, by = .(subjectId, rowname)]
      
    }
    
    # Define a function to run the model for each rowname
    run_model_for_rowname <- function(rowname_value, model) {
      # Filter data for the current rowname
      df_row <- df[rowname == rowname_value]
      
      # Run your model
      model_result <- lmer(model, data = df_row, REML = FALSE)
      
      # Extract summary statistics
      broom.mixed::tidy(model_result) |>
        as.data.frame() |>
        filter(effect == 'fixed' & term != "(Intercept)" & grepl("visitId", term)) |>
        select(term, p.value) |>
        mutate(x = rowname_value) |>
        relocate(x)
    }
    
    # Apply the model to each unique rowname using parallel processing
    results <- mclapply(unique(df$rowname), function(i) {
      tryCatch(run_model_for_rowname(i, model), error = function(e) return(e))
    })
    
    # Remove any with error
    results <- results |>
      purrr::keep( ~ !inherits(.x, 'error')) 
    
    ## Bind
    results <- dplyr::bind_rows(results) |> 
      tidyr::pivot_wider(id_cols = x,
                         names_from = term,
                         values_from = p.value)
    
    return(results)
  }
  
  # Use progressr to monitor progress
  results <- progressr::with_progress({
    p <- progressr::progressor(along = 1:n)
    
    future.apply::future_lapply(1:n, function(x) {
      p()  # Update the progress bar
      shuffle_and_model_iteration(df, model, variables)
    }, future.seed = TRUE)
  })
  
  # Combine all results into a single data frame
  results_df <- bind_rows(results, .id = "replicate_id")
  
  return(results_df)
}

# Example usage:
# results <- permutationTesting(df, n = 500, model = 'value ~ visitId*compliance + age_at_consent + bmi_at_consent + (1 | subjectId)')
# View(results)
