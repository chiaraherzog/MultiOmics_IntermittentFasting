#' FDR Correction within Defined Groups
#'
#' This function applies FDR (False Discovery Rate) correction to p-values within 
#' predefined groups of assays. The groups are defined based on a mapping of 
#' specific assay names to overarching categories.
#'
#' @param df A data frame containing the data for which FDR correction should be applied. 
#' The data frame must include a column `x` that contains the assay names, and one or more 
#' columns containing p-values to be corrected.
#'
#' @return A data frame with FDR-corrected p-values for each group. New columns are added 
#' with the suffix `_adj` to indicate the adjusted p-values.
#'
#' @examples
#' \dontrun{
#' # Assuming `df` is your data frame:
#' result <- FDRcorr(df)
#' }
#'
#' @export
FDRcorr <- function(df, append = T) {
  
  # Define groups of assays for FDR correction
  # Step 1: Define the mapping of assays to overarching groups
  assay_to_group <- c(
    "Blood haemogram" = "clinical",
    "Body composition" = "clinical",
    "Functional sports exam" = "clinical",
    "Vascular and body sonography" = "clinical",
    "Skin histology and transepidermal water loss assay" = "clinical",
    
    "Composite methylation scores: blood" = "blood methylation",
    "Composite methylation scores: buccal" = "buccal methylation",
    "Composite methylation scores: cervical" = "cervical methylation",
    
    "Flow cytometry: T cell staining" = "immune",
    "Flow cytometry: white blood cell staining" = "immune",
    "Immune age: general" = "immune",
    
    "Saliva microbiome: ASVs" = "saliva microbiome asv",
    "Saliva microbiome: ASVs: clr" = "saliva microbiome asv clr",
    "Saliva microbiome: families" = "saliva microbiome family",
    "Saliva microbiome: families: clr" = "saliva microbiome family clr",
    "Saliva microbiome: KO" = "saliva microbiome ko",
    "Saliva microbiome: KO: clr" = "saliva microbiome ko clr",
    
    "Stool microbiome: ASVs" = "stool microbiome asv",
    "Stool microbiome: ASVs: clr" = "stool microbiome asv clr",
    "Stool microbiome: families" = "stool microbiome family",
    "Stool microbiome: families: clr" = "stool microbiome family clr",
    "Stool microbiome: KO" = "stool microbiome ko",
    "Stool microbiome: KO: clr" = "stool microbiome ko clr",
    
    "Saliva nuclear magnetic resonance: normalized" = "saliva metabolome",
    "Urine nuclear magnetic resonance: normalized" = "urine metabolome",
    "Saliva nuclear magnetic resonance: clr" = "urine metabolome clr",
    "Urine nuclear magnetic resonance: clr" = "urine metabolome clr"
  )
  
  # rename variable to x if not yet present
  if("rowname" %in% colnames(df) & !"assay" %in% colnames(df)){
    df <- df |> dplyr::rename(x = rowname)
  }
  
  # Step 2: Apply the mapping to the dataframe
  if('assay' %in% colnames(df)){
    tmp <- df |>
      # Fix for CLR naming
      dplyr::mutate(assay = gsub("_clr_", ": clr_", assay)) |>
      dplyr::mutate(group = assay_to_group[assay])
    
  } else {
    
  tmp <- df |>
    # Fix for CLR naming
    dplyr::mutate(x = gsub("_clr_", ": clr_", x)) |> 
    
    tidyr::separate(x, sep = "_", 
                    into = c('assay', 'x'), 
                    extra = 'merge') |> 
    dplyr::mutate(group = assay_to_group[assay])
  }
  
  # Step 3: Perform FDR correction within groups
  if(append == T){
  tmp <- tmp |>
    dplyr::group_by(group) |> 
    dplyr::mutate(across(contains("p.value"), 
                         ~ p.adjust(., method = 'fdr', n = n()), 
                         .names = "{.col}_adj")) |> 
    dplyr::mutate(across(starts_with("p") & !contains("p.value"), 
                         ~ p.adjust(., method = 'fdr', n = n()), 
                         .names = "{.col}_adj")) |> 
    dplyr::ungroup()
  } else {
    
    # fix if x or rowname
    if("x" %in% colnames(tmp)){
      tmp$rowname <- tmp$x
    }
    
    tmp <- tmp |>
      dplyr::group_by(group) |> 
      dplyr::mutate(across(contains("p.value"), 
                           ~ p.adjust(., method = 'fdr', n = n()))) |> 
      dplyr::mutate(across(starts_with("p") & !contains("p.value"), 
                           ~ p.adjust(., method = 'fdr', n = n()))) |> 
      dplyr::ungroup() |> 
      dplyr::mutate(x = paste0(assay, "_", rowname)) |> dplyr::select(-group)
  }
  
  return(tmp)
}
