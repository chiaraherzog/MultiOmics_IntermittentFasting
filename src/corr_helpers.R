# helper to plot correlations
# plot only strong significant correlations
# output list with all tested correlations
# Author: Charlotte Vavourakis

plot_corr <- function(dat_wide){
  
  out <- list()
  
  corr_all = cor(dat_wide, method = "kendall")
  pval_all = cor.mtest(mat = dat_wide, method = "kendall")$p 
  
  # Filter features based on significance and absolute correlation value
  significant_features = rownames(corr_all)[
    apply((pval_all < 0.01) & (abs(corr_all) >= 0.7) & (pval_all > 0) & (abs(corr_all) < 1), 1, any)
  ]
  corr_all_filtered = corr_all[significant_features, significant_features]
  pval_all_filtered = pval_all[significant_features, significant_features]
  
  # Only color significant correlations
  corrplot(corr_all_filtered,p.mat = pval_all_filtered,
           method = 'circle', type = 'lower',
           sig.level = 0.05, insig='blank',
           tl.col = 'black',tl.srt = 75,
           order = 'AOE', diag = FALSE)$corrPos 
  
  # corrplot(corr_all_filtered,p.mat = pval_all_filtered,
  #          method = 'circle', type = 'lower',
  #          sig.level = 0.01, insig='blank',
  #          tl.col = 'black',
  #          order = 'AOE', diag = FALSE)$corrPos -> p1
  # text(p1$x, p1$y, round(p1$corr, 2))
  
  out[[1]] <- corr_all
  out[[2]] <- pval_all
  
  return(out)
  
}

# helper for filtering and plotting cross-ome feature pairs from precalculated correlation matrices
plot_corr_cross <- function(out){
  
  cross_mat <- function(mat){
    filtered_rows <- grep("salivanmr|urinenmr", rownames(mat), value = TRUE)
    filtered_cols <- grep("stool16Sfam|saliva16Sfam", colnames(mat), value = TRUE)
    mat_filtered <- mat[filtered_rows,filtered_cols, drop = FALSE]
    
    return(mat_filtered)
  }
  
  corr_cross = cross_mat(out[[1]])
  pval_cross = cross_mat(out[[2]])
  
  # Filter features based on significance and absolute correlation value
  significant_rows = rownames(corr_cross)[
    apply((pval_cross < 0.05) & (abs(corr_cross) >= 0.3) & (pval_cross > 0) & (abs(corr_cross) < 1), 1, any)
  ]
  significant_cols = colnames(corr_cross)[
    apply((pval_cross < 0.05) & (abs(corr_cross) >= 0.3) & (pval_cross > 0) & (abs(corr_cross) < 1), 2, any)
  ]
  corr_cross_filtered = corr_cross[significant_rows, significant_cols, drop=FALSE]
  pval_cross_filtered = pval_cross[significant_rows, significant_cols, drop=FALSE]
  
  corrplot(corr_cross_filtered,p.mat = pval_cross_filtered,
           sig.level = 0.05, insig='blank',
           tl.col = 'black',tl.srt = 75)
}

# helper for filtering and plotting metabolome cross-sample feature pairs from precalculated correlation matrices
plot_corr_cross_metabolome <- function(out){
  
  cross_mat <- function(mat){
    filtered_rows <- grep("salivanmr", rownames(mat), value = TRUE)
    filtered_cols <- grep("urinenmr", colnames(mat), value = TRUE)
    mat_filtered <- mat[filtered_rows,filtered_cols, drop = FALSE]
    
    return(mat_filtered)
  }
  
  corr_cross = cross_mat(out[[1]])
  pval_cross = cross_mat(out[[2]])
  
  # Filter features based on significance and absolute correlation value
  significant_rows = rownames(corr_cross)[
    apply((pval_cross < 0.05) & (abs(corr_cross) >= 0.3) & (pval_cross > 0) & (abs(corr_cross) < 1), 1, any)
  ]
  significant_cols = colnames(corr_cross)[
    apply((pval_cross < 0.05) & (abs(corr_cross) >= 0.3) & (pval_cross > 0) & (abs(corr_cross) < 1), 2, any)
  ]
  corr_cross_filtered = corr_cross[significant_rows, significant_cols, drop=FALSE]
  pval_cross_filtered = pval_cross[significant_rows, significant_cols, drop=FALSE]
  
  corrplot(corr_cross_filtered,p.mat = pval_cross_filtered,
           sig.level = 0.05, insig='blank',
           tl.col = 'black',tl.srt = 75)
}

# helper for filtering and plotting microbiome cross-sample feature pairs from precalculated correlation matrices
plot_corr_cross_microbiome <- function(out){
  
  cross_mat <- function(mat){
    filtered_rows <- grep("saliva16Sfam", rownames(mat), value = TRUE)
    filtered_cols <- grep("stool16Sfam", colnames(mat), value = TRUE)
    mat_filtered <- mat[filtered_rows,filtered_cols, drop = FALSE]
    
    return(mat_filtered)
  }
  
  corr_cross = cross_mat(out[[1]])
  pval_cross = cross_mat(out[[2]])
  
  # Filter features based on significance and absolute correlation value
  significant_rows = rownames(corr_cross)[
    apply((pval_cross < 0.05) & (abs(corr_cross) >= 0.3) & (pval_cross > 0) & (abs(corr_cross) < 1), 1, any)
  ]
  significant_cols = colnames(corr_cross)[
    apply((pval_cross < 0.05) & (abs(corr_cross) >= 0.3) & (pval_cross > 0) & (abs(corr_cross) < 1), 2, any)
  ]
  corr_cross_filtered = corr_cross[significant_rows, significant_cols, drop=FALSE]
  pval_cross_filtered = pval_cross[significant_rows, significant_cols, drop=FALSE]
  
  corrplot(corr_cross_filtered,p.mat = pval_cross_filtered,
           sig.level = 0.05, insig='blank',
           tl.col = 'black',tl.srt = 75)
}
