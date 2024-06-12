#' @name renameVarsLME
#' @description
#' Relabel wilcoxon tests for output printing (variables in more legible format)
#' @param out_lmm List output (Wilcoxon tests)

renameVarsLME <- function(
    out_lmm
){
  
  load(here("src/vars.Rdata"))
  
  vars <- vars |> 
    dplyr::mutate(label = ifelse(!is.na(`second name`), `second name`, label))
  
  # Filter and rename vars
  n <- names(out_lmm)
  
  out_lmm <- lapply(n, function(x){
    
    out_lmm[[x]] |> 
      tidyr::separate('x', "_", 
                      into = c("assay", "x"),
                      extra = "merge") |> 
      
      dplyr::select(-starts_with("std.error_")) |> 
      
      dplyr::left_join(dplyr::select(vars, assay, assay2, label, x), by = c('x' = 'x',
                                                                            'assay' = 'assay')) |> 
      dplyr::mutate(variable = ifelse(!is.na(label), label, x),
                    assay = ifelse(!is.na(assay2), assay2, assay)) |> 
      dplyr::select(-any_of(c("x", "label", "assay2"))) |> 
      dplyr::mutate(across(!contains(c("variable", "assay")), ~ signif(., digits = 4))) |> 
      
      dplyr::arrange(assay) |> 
      dplyr::relocate(assay, variable) |> 
      dplyr::distinct() 
    
  })
  
  names(out_lmm) <- n
  
  return(out_lmm)
}
