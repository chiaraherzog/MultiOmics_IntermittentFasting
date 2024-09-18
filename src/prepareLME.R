prepareLME <- function(){
  
  
  load(here("out/out_lmm_factor.Rdata"))
  source(here('src/FDRcorr.R'))
  source(here('src/appendPermutationPval.R'))
  
  load(here('src/vars.Rdata'))
  vars <- vars |>
    dplyr::filter(!grepl("ASVs$|families_clr", assay)) |> 
    dplyr::mutate(x = paste0(assay, "_", x))
  
  
  # Subset relevant models 
  keep <- c("Minimal model", "Basic model with BMI", "Intervention (higher compliance)", "Menopause general", "Menopause interaction")
  out_lmm <- out_lmm[keep]
  
  # FDR Correction
  cat("FDR correction\n")
  out_lmm <- lapply(out_lmm, FDRcorr)
  
  # Unite variables again for next step
  out_lmm <- lapply(out_lmm, function(i){
    i |> tidyr::unite(col = "x", assay:x, sep = "_")
  })
  
  cat("Append permutation: ")
  load(here("out/outPermMinimal.Rdata"))
  load("out/outPerm.Rdata")
  load(here("out/outPermIntervention.Rdata"))
  
  cat("minimal ")
  out_lmm$`Minimal model` <- appendPermutationComparisonPval_OutPut(observed = out_lmm$`Minimal model`,
                                                                    random = outPermMinimal,
                                                                    variables = vars$x)
  cat("[done]; basic ")
  
  out_lmm$`Basic model with BMI` <- appendPermutationComparisonPval_OutPut(observed = out_lmm$`Basic model with BMI`,
                                                                    random = outPerm,
                                                                    variables = vars$x)
  cat("[done]; intervention ")
  
  out_lmm$`Intervention (higher compliance)` <- appendPermutationComparisonPval_OutPut(observed = out_lmm$`Intervention (higher compliance)`,
                                                                           random = outPermIntervent,
                                                                           variables = vars$x)
  cat("[done];\n")
  
  # Relabel for outputs
  cat("Rename variables for easier understanding\n")
  source(here('src/renameVarsLME.R'))
  out_lmm <- renameVarsLME(out_lmm)
  
  
  # Tidy up
  cat("Additional tidying\n")
  out_lmm <- lapply(out_lmm, function(i){
    i |> dplyr::relocate(assay, group, variable) |> 
      dplyr::rename_at(vars(contains("_adj")), ~ gsub("_adj", "_fdr.corr", .)) |> 
      dplyr::filter(!group %in% c("urine metabolome clr",
                                  "saliva metabolome clr",
                                  "stool microbiome asv",
                                  "saliva microbiome asv",
                                  "saliva microbiome family clr",
                                  "stool microbiome family clr"))
  })
  
  
  
  return(out_lmm)

  
}