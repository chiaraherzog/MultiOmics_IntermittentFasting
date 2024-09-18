#' Get permutation cmoparison p value on multiple columns
#'
#' This function calculates empirical p-values for multiple variables and columns
#' using values from permutation testing. It applies the permutation test across
#' all relevant columns in the `observed` data frame that match the pattern "p.value"
#' and "visitId", comparing them to corresponding columns in the `random` data frame.
#'
#' @param observed A data frame containing the observed p-values for each variable.
#' @param random A data frame containing the permutation results for each variable.
#' @param variables A vector of variables for which to calculate empirical p-values.
#'
#' @return A data frame in wide format, where each row corresponds to a variable,
#' and each column corresponds to an empirical p-value calculated for the specified columns.
#' @examples
#' # Assume `observed` and `random` are data frames with appropriate structure
#' # variables <- unique(observed$x)
#' # results <- perform_permutation_testing(observed, random, variables)
#' # print(results)
#' @export

getPermutationComparisonPval <- function(observed, random, variables) {
  
  # Identify columns to permute (assuming "p.value" and "visitId" columns are relevant)
  columnsToPermute <- colnames(observed)[grepl("p.value", colnames(observed)) & grepl("visitId", colnames(observed))]
  
  # Define the permutation testing function
  permutationResult <- function(observed, random, variable, column){
    
    # Extract the observed p-value for the specific variable and column
    pval_obs <- observed[observed$x == variable, ][[column]]
    
    # Extract the permutation results for the specific variable and column
    perm <- random[random$x == variable, ][[gsub("p.value_", "", column)]]
    
    # Calculate the empirical p-value
    p <- sum(abs(perm) <= pval_obs) / length(perm)
    
    return(p)
  }
  
  # Create a data frame of all combinations of variables and columns
  combinations <- tidyr::expand_grid(variable = variables, column = columnsToPermute)
  
  # Apply the permutationResult function to all combinations
  results <- combinations |> 
    mutate(empirical_p = purrr::map2_dbl(variable, column, ~ permutationResult(observed, random, .x, .y)))
  
  # Pivot the results to wide format
  results_wide <- results |> tidyr::pivot_wider(names_from = column, values_from = empirical_p) |> 
    dplyr::rename(x = variable)
  
  # Now, we check: if p value is below 0.05 in results wide -> significant and we take forward the significant p value from the original frame; if not, p value is set to pvalue from results_wide (as a dummy)
  
  permuted <- observed |>
    dplyr::left_join(results_wide, by = 'x') |> 
    rowwise() |>
    dplyr::mutate(across(ends_with(".x"),
                         ~ ifelse(get(cur_column()) < 0.05, .x, 0.5),
                         .names = "{col}"))
  
  # Clean up column names (remove the ".x" suffix)
  permuted_clean <- permuted |> 
    dplyr::select(-ends_with(".y")) |> 
    dplyr::rename_with(~ gsub("\\.x$", "", .), ends_with(".x"))
  
  return(permuted_clean)
}
