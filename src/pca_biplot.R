pca_biplot <- function(pca.plot, pcObject, components=1:2, scale=1, obs.scale = 1 - scale,textsize=3, contributions=NULL, arrow.colour="#1b69a1",text.colour='black', label.colour=NULL, arrow.alpha=0.5, label.alpha=0.8, varname.adjust=1.25,var.length=3){
  #' PCA Biplot Function - B. Theeuwes for TG
  #'
  #' Generates a biplot for Principal Component Analysis (PCA) results.
  #'
  #' @param pca.plot The PCA plot made with ggplot
  #' @param pcObject The principal components object, i.e. output from prcomp.
  #' @param components The components to plot (default: 1:2).
  #' @param scale The scale factor for variables (default: 1).
  #' @param obs.scale The scale factor for observations (default: 1 - scale).
  #' @param textsize The size of text labels (default: 3).
  #' @param contributions Number of contributions, that is, number of arrows plotted. (default: NULL (all)).
  #' @param arrow.colour The color of arrows (default: "#1b69a1").
  #' @param text.colour The textcolor of text labels (default: 'black').
  #' @param label.colour The background color of text labels (default: NULL, means no label is plotted).
  #' @param arrow.alpha The transparency of arrows (default: 0.5).
  #' @param label.alpha The transparency of labels (default: 0.8).
  #' @param varname.adjust Adjustment factor the placement of the variable names, >= 1 means farther from the arrow (default: 1.25).
  #' @param var.length Length of the arrows (default: 3).
  #'
  #' @return A ggplot2 object representing the biplot.
  #'
  #' @details This function adds arrows representing the loadings of the original variables to an existing ggplot. The function allows customization of several parameters such as the scaling of the data, size and colors of the arrows, and text properties.
  #'
  #' @examples
  #' # Example usage:
  #' pca_biplot(pca.plot = figure, pcObject = pcs)
  #'
  
  library(ggbiplot)
  # Extracting labels from principal components object
  labels=rownames(pcObject$rotation) 
  
  if (length(components) > 2) {
    warning("components = ", components, " is not of length 2. Only the first 2 will be used")
    components <- components[1:2]
  }
  
  # Calculating singular value decomposition
  svd <- get_SVD(pcObject)
  n <- svd$n
  d <- svd$D
  u <- svd$U
  v <- svd$V
  
  # Scaling
  nobs.factor <- ifelse(inherits(pcObject, "prcomp"), sqrt(n - 
                                                          1), sqrt(n))
  angle <- circle_chol <- ed <- hjust <- mu <- sigma <- varname <- xvar <- yvar <- NULL
  components <- pmin(components, ncol(u))
  df.u <- as.data.frame(sweep(u[, components], 2, d[components]^obs.scale, 
                              FUN = "*"))
  v <- sweep(v, 2, d^1, FUN = "*")
  df.v <- as.data.frame(v[, components])
  names(df.u) <- c("xvar", "yvar")
  names(df.v) <- names(df.u)
  
  # Scaling factor
  df.u <- df.u * nobs.factor
  
  # Scaling for plotting
  r <- sqrt(qchisq(0.68, df = 2)) * prod(colMeans(df.u^2))^(1/4)
  v.scale <- rowSums(v^2)
  df.v <- r * df.v/sqrt(max(v.scale))
  
  # Calculating angles and adjustments for plotting
  df.v$angle <- with(df.v, (180/pi) * atan(yvar/xvar))
  df.v$hjust = with(df.v, (1 - varname.adjust * sign(xvar))/2)
  
  # If number of contributions to plot is given, subset the df
  if (length(contributions)>0){
    df.v <- df.v[1:contributions,]
    }
  
  
  # Adding arrows to the plot
  arrow_style <- arrow(length = unit(1/3, "picas"), type = "closed", 
                       angle = 15)
  pca.plot <- pca.plot + geom_segment(data = df.v, aes(x = 0, y = 0, 
                                                       xend = xvar*var.length, yend = yvar*var.length), arrow = arrow_style, 
                                      color = arrow.colour, linewidth = 1,alpha = arrow.alpha)
  # Adding labels to the arrows
  
  if (length(label.colour>0)){
  
  pca.plot <- pca.plot + 
    # geom_text(data = df.v, aes(label = rownames(df.v), 
    #           x = xvar*var.length, y = yvar*var.length, hjust = hjust), 
    #           color = text.colour, size = textsize) +
    geom_label(alpha = label.alpha, color = text.colour,data = df.v, aes(label = rownames(df.v), 
                                                                       x = xvar*var.length, y = yvar*var.length, hjust = hjust), 
               color = text.colour, size = textsize)
  }else{
    pca.plot <- pca.plot + 
      geom_text(data = df.v, aes(label = rownames(df.v),
                x = xvar*var.length, y = yvar*var.length, hjust = hjust),
                color = text.colour, size = textsize) 

    
    
  }
  
  
  return(pca.plot)
}