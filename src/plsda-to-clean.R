load("data/data_normalized.Rdata")
vars <- read.table("src/clinical_variables.csv", header = T, sep = ',') |> dplyr::pull(x)
labels = read.table("src/clinical_variables.csv", header = T, sep = ',')
labels <-  labels |> 
  dplyr::add_row(x = 'age_at_consent',
                 label = 'age at consent')

data_pca <- as.data.frame(wideFormat(data[vars,,],
                                     colData = c('subjectId', 'interventionId', 'time', 'compliance', 'age_at_consent'))) |> 
  dplyr::filter(interventionId != "S" & compliance == 'high' & !time %in% c(12, 18) & time %in% c(0, 6)) |> 
  tibble::column_to_rownames('primary')

# remove superfluous parts of rownames and relabel
colnames(data_pca) <- gsub("Functional.sports.exam_|Vascular.and.body.sonography_|Body.composition_|Blood.haemogram_|Skin.histology.and.transepidermal.water.loss.assay_", "", colnames(data_pca))

data_pca <- data_pca |> 
  dplyr::rename_at(vars(labels[labels$x %in% colnames(data_pca),]$x), ~ labels[labels$x %in% colnames(data_pca),]$label)


library(mixOmics)

X <- data_pca |> dplyr::select(-c(subjectId, interventionId, time, compliance, 'age at consent'))
Y <- data_pca$time
result.plsda <- plsda(X, Y,ncomp = 5) # run the method
plotIndiv(result.plsda) # plot the samples
plotVar(result.plsda,overlap = F) # plot the variables

plotLoadings(result.plsda, comp = 1, method = 'mean', contrib = 'max')
plotLoadings(result.plsda, comp = 2, method = 'mean', contrib = 'max')
plotLoadings(result.plsda, comp = 3, method = 'mean', contrib = 'max')
plotLoadings(result.plsda, comp = 4, method = 'mean', contrib = 'max')
plotLoadings(result.plsda, comp = 5, method = 'mean', contrib = 'max')

result.plsda$loadings



library(mixOmics)
library(ggplot2)
library(reshape2)

# Function to perform PLS-DA, plot biplot, and plot loadings
plsda_analysis <- function(X, Y, ncomp = 5, top_n = 10) {
  
  # Run the PLS-DA
  result.plsda <- plsda(X, Y, ncomp = ncomp)
  
  # plotting individuals
  plot_indiv <- function(result, comp1 = 1, comp2 = 2, design) {
    # Extract scores for the specified components
    scores <- as.data.frame(result$variates$X[, c(comp1, comp2)])
    scores$DesignFactor <- design
    
    # Extract loadings for the specified components
    loadings <- as.data.frame(result$loadings$X[, c(comp1, comp2)])
    loadings$Variable <- rownames(loadings)
    
    # Create the ggplot
    p <- ggplot(scores, aes(x = scores[,1], y = scores[,2], color = as.factor(DesignFactor))) +
      geom_point(size = 2, alpha = 0.8) +
      labs(x = paste("Component", comp1), y = paste("Component", comp2),
           title = paste("PLS-DA Plot of Components", comp1, "and", comp2)) +
      theme_bw() +
      theme(legend.title = element_blank()) +
      stat_ellipse() +
      scale_colour_manual(values = cols[c(7, 3)])
    
    return(p)
  }
  
  plotIndiv <- plot_indiv(result.plsda, design = Y)
    
  # plot loadings
  plot_loadings <- function(result.plsda, comp = 1, Y, method = 'mean', contrib = 'max', top_n = top_n) {
    # Extract the loadings for the specified component
    loadings <- as.data.frame(result$loadings$X[, comp])
    loadings$Variable <- rownames(loadings)
    colnames(loadings) <- c("Loading", "Variable")
    
    # Create a data frame to store the calculated contributions for each variable
    contributions <- data.frame(Variable = loadings$Variable, GroupContrib = NA)
    
    # For each variable, calculate the mean (or median) loading per group in Y
    for (var in loadings$Variable) {
      # Get the values of the variable across samples
      var_values <- result$X[, var]
      
      # Calculate the mean (or method) per group
      group_means <- aggregate(var_values, by = list(Y), FUN = method, na.rm = T)
      
      # Find the group with the maximum (or minimum) mean
      if (contrib == 'max') {
        max_group <- group_means[which.max(group_means$x), 1]
      } else {
        max_group <- group_means[which.min(group_means$x), 1]
      }
      
      # Store the group that has the maximum contribution
      contributions$GroupContrib[contributions$Variable == var] <- max_group
    }
    
    # Merge the loadings and contributions data frames
    loadings <- merge(loadings, contributions, by = "Variable")
    
    # Select the top_n variables based on absolute loading values
    top_loadings <- loadings[order(-abs(loadings$Loading)), ][1:top_n, ]
    
    # Create a barplot using ggplot2
    ggplot(top_loadings, aes(x = reorder(Variable, abs(Loading)), y = Loading, fill = as.factor(GroupContrib))) +
      geom_bar(stat = "identity") +
      coord_flip() +
      labs(x = "", y = "Loadings", 
           title = paste("Component", comp)) +
      theme_bw() +
      theme(panel.grid = element_blank(),
            axis.ticks.y = element_blank()) +
      scale_fill_manual(values =cols[c(3, 7)],
                        name = '')
  }
  
  loadings1 <- plot_loadings(plot_loadings, comp = 1, Y = Y, top_n = top_n)
  loadings2 <- plot_loadings(plot_loadings, comp = 2, Y = Y, top_n = top_n)
  loadings3 <- plot_loadings(plot_loadings, comp = 3, Y = Y, top_n = top_n)
  loadings4 <- plot_loadings(plot_loadings, comp = 4, Y = Y, top_n = top_n)
  
  return(list(result = result.plsda,
              plotIndiv = plotIndiv,
              loadings1 = loadings1,
              loadings2 = loadings2,
              loadings3 = loadings3,
              loadings4 = loadings4))
  
}

out <- plsda_analysis(X, Y)


out$plotIndiv
out$loadings2

p2 <- (out$loadings1/out$loadings2/out$loadings3/out$loadings4) + plot_layout(guides = 'collect')

library(patchwork)
(out$plotIndiv | p2) + plot_layout(widths = c(1.5, 0.5))

