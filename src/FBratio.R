
# calc_FBratio : calculate the Firmicutes:Bacteroidota ratio for an ASV/OTU-table (wide format)
# plot_FBratio : boxplot changes over time, test significance
# calc_CBratio : calculate the Clostridia:Bacteroidia ratio for an ASV/OTU-table (wide format)
# plot_CBratio : boxplot changes over time, test significance
# Author: Charlotte Vavourakis

# helper function for outlier removal function
remove_outliers <- function(data, variable) {
  q1 <- quantile(data[[variable]], 0.25)
  q3 <- quantile(data[[variable]], 0.75)
  iqr <- q3 - q1
  lower_bound <- q1 - 1.5 * iqr
  upper_bound <- q3 + 1.5 * iqr
  data <- data %>%
    filter(get(variable) >= unname(lower_bound), get(variable) <= unname(upper_bound))
  return(data)
}

calc_FBratio <- function(ASVdat){
  
  # make a long format, FB only
  datFB <- ASVdat %>%
    filter(Phylum %in% c("Firmicutes","Bacteroidota")) %>%
    select(-Kingdom,-Class,-Order,-Family,-Genus,-Species) %>%
    pivot_longer(cols = contains("FM"), names_to = "sampleId", values_to = "count")
  
  # sum counts for each phylum
  datF <- datFB %>%
    filter(Phylum == "Firmicutes") %>%
    group_by(sampleId) %>%
    summarise(total_firmicutes_counts = sum(count))
  datB <- datFB %>%
    filter(Phylum == "Bacteroidota") %>%
    group_by(sampleId) %>%
    summarise(total_bacteroidota_counts = sum(count))
  
  # calculate the ratio
  datFB <- inner_join(datF, datB, by = "sampleId") %>%
    mutate(FB_ratio = total_firmicutes_counts / total_bacteroidota_counts)
  
  return(datFB)
  
}

plot_FBratio <- function(dat, remove.outliers = F, subset) {
  
  # make long format
  df = dat %>%
    dplyr::select(subjectId,visitId, FB_ratio) %>%
    na.omit()
  
  if (remove.outliers == T){
    df <- remove_outliers(df, "FB_ratio")
  }
  
  # keep only participants with complete time series for significance testing
  # Count the number of unique visitIds for each subjectId
  visit_counts <- df %>%
    group_by(subjectId) %>%
    summarise(n_visits = n_distinct(visitId))
  
  # Filter subjectIds with all four visitIds
  complete_subjects <- visit_counts %>%
    filter(n_visits == 4) %>%
    pull(subjectId)
  
  # Filter the dataframe to keep only the complete subjectIds
  df <- df %>%
    filter(subjectId %in% complete_subjects)
  
  compare = list( c("M0", "M2"), c("M0", "M4"), c("M0", "M6") )
  
  df$visitId = factor(df$visitId, levels=c("M0","M2","M4","M6"))
  
  plot = ggboxplot(df, x = "visitId", y = "FB_ratio") + 
    geom_line(aes(group = subjectId), colour = "#1b69a1", alpha = 0.3) + 
    stat_compare_means(comparisons = compare, paired=TRUE)+ 
    stat_compare_means(label.y = max(df$FB_ratio)-0.5,
                       label.x = 1.8) +
    ylab("Firmicutes:Bacteroidota ratio") +
    xlab("") +
    theme(legend.position="none") +
    ggtitle(subset)
  
  return(plot)
}

calc_CBratio <- function(ASVdat){
  
  # make a long format, CB only
  datCB <- ASVdat %>%
    filter(Class %in% c("Clostridia","Bacteroidia")) %>%
    select(-Kingdom,-Phylum,-Order,-Family,-Genus,-Species) %>%
    pivot_longer(cols = contains("FM"), names_to = "sampleId", values_to = "count")
  
  # sum counts for each phylum
  datC <- datCB %>%
    filter(Class == "Clostridia") %>%
    group_by(sampleId) %>%
    summarise(total_clostridia_counts = sum(count))
  datB <- datCB %>%
    filter(Class ==  "Bacteroidia") %>%
    group_by(sampleId) %>%
    summarise(total_bacteroidia_counts = sum(count))
  
  # calculate the ratio
  datCB <- inner_join(datC, datB, by = "sampleId") %>%
    mutate(CB_ratio = total_clostridia_counts / total_bacteroidia_counts)
  
  return(datCB)
  
}

plot_CBratio <- function(dat, remove.outliers = F, subset) {

  # make long format
  df = dat %>%
    dplyr::select(subjectId,visitId, CB_ratio) %>%
    na.omit()
  
  if (remove.outliers == T){
    df <- remove_outliers(df, "CB_ratio")
  }
  
  # keep only participants with complete time series for significance testing
  # Count the number of unique visitIds for each subjectId
  visit_counts <- df %>%
    group_by(subjectId) %>%
    summarise(n_visits = n_distinct(visitId))
  
  # Filter subjectIds with all four visitIds
  complete_subjects <- visit_counts %>%
    filter(n_visits == 4) %>%
    pull(subjectId)
  
  # Filter the dataframe to keep only the complete subjectIds
  df <- df %>%
    filter(subjectId %in% complete_subjects)
  
  compare = list( c("M0", "M2"), c("M0", "M4"), c("M0", "M6") )
  
  df$visitId = factor(df$visitId, levels=c("M0","M2","M4","M6"))
  
  plot = ggboxplot(df, x = "visitId", y = "CB_ratio") + 
    geom_line(aes(group = subjectId), colour = "#1b69a1", alpha = 0.3) + 
    stat_compare_means(comparisons = compare, paired=TRUE)+ 
    stat_compare_means(label.y = max(df$CB_ratio)-0.5,
                       label.x = 1.8) +
    ylab("Clostridia:Bacteroidia ratio") +
    xlab("") +
    theme(legend.position="none") +
    ggtitle(subset)
  
  return(plot)
}