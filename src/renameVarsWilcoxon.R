#' @name renameVarsWilcoxon
#' @description
#' Relabel wilcoxon tests for output printing (variables in more legible format)
#' @param out_wilcoxon List output (Wilcoxon tests)

renameVarsWilcoxon <- function(
    out_wilcoxon
){
  
  load(here("src/vars.Rdata"))
  
  vars <- vars |> 
    dplyr::mutate(label = ifelse(!is.na(`second name`), `second name`, label)) |> 
    dplyr::select(x, label, assay, assay2)
  
  n <- names(out_wilcoxon)
  
  out_wilcoxon <- lapply(n, function(x){
    out_wilcoxon[[x]] |>
      tidyr::separate(rowname, into = c('assay', "variable"), sep = "_", extra = 'merge') |>
      dplyr::left_join(vars, by = c('variable' = 'x',
                                    'assay' = 'assay')) |>
      dplyr::mutate(variable = ifelse(!is.na(label), label, variable),
                    assay = ifelse(!is.na(assay2), assay2, assay)) |>
      dplyr::select(-any_of(c("x", "label"))) |>
      dplyr::arrange(assay, pick(any_of(c("p_M6", "p value (higher compliance versus other)", "p value (I versus K, high compliance)")))) |>
      dplyr::select(-any_of(c("population name", "name", "fixable viability dye and antibodies α-", "staining", "type", "antibodies", "targets", "main_analysis",
                              "assay2"))) |>
      dplyr::mutate(across(contains(c("p_M2", "p_M4", "p_M6", "p value (high versus other)", "p value (I versus K, high compliance)", "sd_")), ~ signif(., digits = 4)))
  })
  names(out_wilcoxon) <- n
  
  return(out_wilcoxon)
}
