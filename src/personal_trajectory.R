personal_trajectories <- function(data, subjects, feature){
  
  # Initiate outputs -----
  df_personal_r <- matrix(ncol = length(unique(subjects)),
                          nrow = 1)
  df_personal_p <- matrix(ncol = length(unique(subjects)),
                          nrow = 1)
  df_personal_slope <- matrix(ncol = length(unique(subjects)),
                              nrow = 1)
  df_personal_wald <- matrix(ncol = length(unique(subjects)),
                             nrow = 1)
  
  colnames(df_personal_r) <- unique(subjects)
  colnames(df_personal_p) <- unique(subjects)
  colnames(df_personal_slope) <- unique(subjects)
  colnames(df_personal_wald) <- unique(subjects)
  
  rownames(df_personal_r) <- feature
  rownames(df_personal_p) <- feature
  rownames(df_personal_slope) <- feature
  rownames(df_personal_wald) <- feature
  
  for(i in subjects){
    # subset df
    tmp <- data[data$rowname==feature & data$subjectId==i,]
    tmp <- tmp[!is.na(tmp$value),]
    
    if(nrow(tmp) >= 3){
      # correlation
      r <- cor.test(tmp$value, tmp$time,method = 'pearson')
      
      # Rfit
      f <- tryCatch(Rfit::rfit(value ~ time, data = tmp), error = function(error) NULL)
      
      # Overall wald test
      wald <- tryCatch(Rfit::wald.test.overall(f), error = function(error) NULL)
      # wald <- tryCatch(aod::wald.test(Sigma = vcov(f), b = coef(f), Terms = 1:2), error = function(error) NULL)
      
      df_personal_r[,i] <- as.numeric(r$estimate)
      df_personal_p[,i] <- ifelse(!exists("f"), NA, as.numeric(r$p.value))
      df_personal_slope[,i] <- ifelse(!exists("f"), NA, as.numeric(f$coefficients[2]))
      # df_personal_wald[,i] <- ifelse(!exists("wald"), NA, as.numeric(wald$result$chi2[3]))
      df_personal_wald[,i] <- ifelse(!exists("wald"), NA, as.numeric(wald$p.value))
      
    } else if(nrow(tmp)==2 & all(c(0, 6) %in% tmp$time)){
        df_personal_slope[,i] <- (tmp[tmp$time==6,]$value - tmp[tmp$time==0,]$value)/6
    }
  }
  
  # Put together output
  list <- list(r = df_personal_r,
               p = df_personal_p,
               slope = df_personal_slope,
               wald = df_personal_wald)
  
  # list <- lapply(list, filter_nas)
  
  # Put together an output
  return(list)
  
}

filter_nas <- function(output){
  
  tmp <- output |> 
    as.data.frame() |> 
    dplyr::select_if( ~ !all(is.na(.))) |> 
    dplyr::select(-featureid)
  
  return(tmp)
  
}
