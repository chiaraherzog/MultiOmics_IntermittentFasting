#' @name t_tests_sep
#' @description
#' Run t tests (per intervention) - no comparison of I versus K.
#'  description
#' @param dat Long dataframe of change from baseline data (to compute change values)
#' @param dat_raw Long dataframe of raw values to compute e.g. baseline values
#' @param intervention Study arm: I or S. If I, will compute between-intervention results

t_tests_sep <- function(dat, dat_raw){
  
  cat('Formatting data... \n')
  # filter complete cases
  complete_cases <- dat |>
    dplyr::filter(!is.na(value)) |>
    dplyr::select(assay, subjectId, visitId, rowname) |>
    dplyr::distinct() |>
    dplyr::group_by(assay, subjectId, rowname) |>
    dplyr::filter(all(c("M0", "M2", "M4", "M6") %in% visitId)) |>
    dplyr::pull(subjectId) |> unique()
  
  dat <- dat |> dplyr::filter(subjectId %in% complete_cases)
  dat_raw <- dat_raw |> dplyr::filter(subjectId %in% complete_cases)  
  
  # Pivot wide for easier filtering
  dat <- dat |> tidyr::pivot_wider(id_cols = c('assay', 'rowname', 'interventionId', 'subjectId', 'compliance'),
                                   names_from = 'visitId',
                                   values_from = 'value')
  
  dat_raw <- dat_raw |> tidyr::pivot_wider(id_cols = c('assay', 'rowname', 'interventionId', 'subjectId', 'compliance'),
                                           names_from = 'visitId',
                                           values_from = 'value')
  
  l = as.list(unique(dat_raw$rowname))
  
  # -------------------------------------------------
  cat('Compute baseline values (all)...\n')
  
  pb <- txtProgressBar(min = 0, max = length(l), style = 3,
                       width = 50)
  
  baseline_all <- do.call(rbind, (lapply(l, function(x){
    
    setTxtProgressBar(pb, which(l == x))
    if(grepl("sports|histology|sonography", x)){
      dat_raw |> 
        dplyr::filter(rowname == x) |> 
        dplyr::filter(!is.na(M0) & !is.na(M6)) |> 
        dplyr::reframe(rowname = x,
                       mean_baseline = mean(M0),
                       sd_baseline = sd(M0)) |> 
        dplyr::ungroup() 
    } else {
      dat_raw |> 
        dplyr::filter(rowname == x) |> 
        dplyr::filter(!is.na(M0) & !is.na(M2) & !is.na(M4) & !is.na(M6)) |> 
        dplyr::reframe(rowname = x,
                       mean_baseline = mean(M0),
                       sd_baseline = sd(M0)) |> 
        dplyr::ungroup() 
    }
    
  })))
  
  # -------------------------------------------------
  cat('\nCompute baseline values (high compliance)...\n')
  pb <- txtProgressBar(min = 0, max = length(l), style = 3,
                       width = 50)
  
  baseline_high <- do.call(rbind, (lapply(l, function(x){
    
    setTxtProgressBar(pb, which(l == x))
    
    if(grepl("sports|histology|sonography", x)){
      dat_raw |> 
        dplyr::filter(rowname == x & compliance == 'high') |> 
        dplyr::filter(!is.na(M0) & !is.na(M6)) |> 
        dplyr::reframe(rowname = x,
                       mean_baseline = mean(M0),
                       sd_baseline = sd(M0)) |> 
        dplyr::ungroup() 
    } else {
      dat_raw |> 
        dplyr::filter(rowname == x & compliance == 'high') |> 
        dplyr::filter(!is.na(M0) & !is.na(M2) & !is.na(M4) & !is.na(M6)) |> 
        dplyr::reframe(rowname = x,
                       mean_baseline = mean(M0),
                       sd_baseline = sd(M0)) |> 
        dplyr::ungroup() 
    }
    
  })))
  
  # -------------------------------------------------
  cat('\nCompute change / t tests...\n')
  pb <- txtProgressBar(min = 0, max = length(l), style = 3,
                       width = 50)
  
  change_all <- do.call(rbind, lapply(l, function(x){
    
    setTxtProgressBar(pb, which(l == x))
    
    if(grepl("sports|histology|sonography", x)){
      
      tmp <- dat |> 
        dplyr::filter(rowname == x) |> 
        dplyr::filter(!is.na(M0) & !is.na(M6)) |> 
        dplyr::select(c(M0, M2, M4, M6))
      
      p_m6 <- wilcox.test(tmp$M0, tmp$M6, paired = T)$p.value
      
      tmp |> 
        dplyr::reframe(rowname = x,
                       dM2 = NA,
                       sd_M2 = NA,
                       p_M2 = NA,
                       dM4 = NA,
                       sd_M4 = NA,
                       p_M4 = NA,
                       dM6 = mean(M6),
                       sd_M6 = sd(M6),
                       p_M6 = p_m6)
    } else {
      
      tmp <- dat |> 
        dplyr::filter(rowname == x) |> 
        dplyr::filter(!is.na(M0) & !is.na(M2) & !is.na(M4) & !is.na(M6)) |> 
        dplyr::select(c(M0, M2, M4, M6))
      
      p_m2 <- wilcox.test(tmp$M0, tmp$M2, paired = T)$p.value
      p_m4 <- wilcox.test(tmp$M0, tmp$M4, paired = T)$p.value
      p_m6 <- wilcox.test(tmp$M0, tmp$M6, paired = T)$p.value
      
      tmp |> 
        dplyr::reframe(rowname = x,
                       dM2 = mean(M2),
                       sd_M2 = sd(M2),
                       p_M2 = p_m2,
                       dM4 = mean(M4),
                       sd_M4 = sd(M4),
                       p_M4 = p_m4,
                       dM6 = mean(M6),
                       sd_M6 = sd(M6),
                       p_M6 = p_m6)
    }
    
  }))
  
  # -------------------------------------------------
  cat('\nCompute change / t tests in high compliance group only...\n')
  pb <- txtProgressBar(min = 0, max = length(l), style = 3,
                       width = 50)
  
  change_high <- do.call(rbind, lapply(l, function(x){
    
    setTxtProgressBar(pb, which(l == x))
    
    if(grepl("sports|histology|sonography", x)){
      
      tmp <- dat |> 
        dplyr::filter(rowname == x & compliance == 'high') |> 
        dplyr::filter(!is.na(M0) & !is.na(M6)) |> 
        dplyr::select(c(M0, M2, M4, M6))
      
      p_m6 <- wilcox.test(tmp$M0, tmp$M6, paired = T)$p.value
      
      tmp |> 
        dplyr::reframe(rowname = x,
                       dM2 = NA,
                       sd_M2 = NA,
                       p_M2 = NA,
                       dM4 = NA,
                       sd_M4 = NA,
                       p_M4 = NA,
                       dM6 = mean(M6),
                       sd_M6 = sd(M6),
                       p_M6 = p_m6)
      
    } else {
      
      tmp <- dat |> 
        dplyr::filter(rowname == x & compliance == 'high') |> 
        dplyr::filter(!is.na(M0) & !is.na(M2) & !is.na(M4) & !is.na(M6)) |> 
        dplyr::select(c(M0, M2, M4, M6))
      
      p_m2 <- wilcox.test(tmp$M0, tmp$M2, paired = T)$p.value
      p_m4 <- wilcox.test(tmp$M0, tmp$M4, paired = T)$p.value
      p_m6 <- wilcox.test(tmp$M0, tmp$M6, paired = T)$p.value
      
      tmp |> 
        dplyr::reframe(rowname = x,
                       dM2 = mean(M2),
                       sd_M2 = sd(M2),
                       p_M2 = p_m2,
                       dM4 = mean(M4),
                       sd_M4 = sd(M4),
                       p_M4 = p_m4,
                       dM6 = mean(M6),
                       sd_M6 = sd(M6),
                       p_M6 = p_m6)
    }
  }))
  
  cat('\nReformatting...\n')
  
  colnames(change_all) <- gsub("^d", "", colnames(change_all))
  overall <- baseline_all |> 
    dplyr::left_join(change_all, by = 'rowname')
  
  colnames(change_high) <- gsub("^d", "", colnames(change_high))
  high <- baseline_high |> 
    dplyr::left_join(change_high, by = 'rowname')
  
  cat('Compiling...\n')
  out <- list(overall = overall,
              `high compliance only` = high)
  
  cat('Returning\n')
  return(out)
}
