#' @name renameVarsrmcorr
#' @description
#' Relabel RMCorr for output printing (variables in more legible format)
#' @param rmcorr List output (RMcorr tests)

renameVarsrmcorr <- function(
    rmcorr
){
  
  load(here("src/vars.Rdata"))
  
  vars <- vars |> 
    dplyr::mutate(label = ifelse(!is.na(`second name`), `second name`, label)) |> 
    dplyr::select(label, assay, assay2, x)
  
  rmcorr <- rmcorr |> 
    dplyr::select(-assay2) |> 
    tidyr::separate(measure1, "_", into = c("assay", "x"),
                    extra = 'merge') |> 
    dplyr::mutate(assay = gsub("-log$", "", assay)) |> 
    dplyr::left_join(vars, by = c('x' = 'x',
                                  'assay' = 'assay')) |>
    dplyr::mutate(variable1 = ifelse(!is.na(label), label, x),
                  assay1 = ifelse(!is.na(assay2), assay2, assay)) |>
    dplyr::select(-c(assay2, label, assay, x)) |> 
    tidyr::separate(measure2, "_", into = c("assay", "x"),
                    extra = 'merge') |> 
    dplyr::mutate(assay = gsub("-log$", "", assay)) |> 
    dplyr::left_join(vars, by = c('x' = 'x',
                                  'assay' = 'assay')) |>
    dplyr::mutate(variable2 = ifelse(!is.na(label), label, x),
                  assay2 = ifelse(!is.na(assay2), assay2, assay)) |>
    dplyr::select(-c(label, assay, x)) |> 
    dplyr::relocate(assay1, variable1, assay2, variable2)
  
  return(rmcorr)
}
