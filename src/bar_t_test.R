#' Make bar plot of t-test analysis on metabolite data
#'
#' This function performs analysis on metabolite data and creates a bar plot showing the proportion of significant metabolites at different time points.
#'
#' @param t_results A data frame containing t-test results with columns 'variable', 'p_M2', 'p_M4', and 'p_M6'.
#' @param sampletype A character specifying the type of sample. Default is "urine". Options are "urine" or "saliva".
#' @param class A character specifying the class of metabolites to analyze. Default is "Sub.class". Options are "Super.class", "Main.class", or "Sub.class".
#' @param plot.reference Logical indicating whether to plot the general division of metabolites on the left. Default is TRUE.
#' @param legend.name Name of the legend. Default is "Sub class".
#'
#' @return A ggplot object representing the bar plot showing the proportion of significant metabolites at different time points.
#' 
#' @examples
#' # Example usage:
#' bar_t_test(t_results, sampletype = "urine", class = "Sub.class", plot.reference = TRUE)
#'
#' @import ggplot2
#' @import tidyr
#' @import dplyr
#' @importFrom stats readRDS
#' @export
bar_t_test <- function(t_results, sampletype="urine", class="Sub.class", plot.reference=T, legend.name="Sub class"){

  

# Load the correct reference map  
if (sampletype == "urine") {
RefMet_mapped <-  readRDS("out/RefMet_mapped_urine.Rds")
} else if (sampletype == "saliva") {
RefMet_mapped <-  readRDS("out/RefMet_mapped_saliva.Rds")
}
  
cols <- c("Nicotinic acid alkaloids" = "#1b69a1", "Amino acids and peptides" = "#48a0af", "-" = "#71b5a9", 
              "Tryptophan alkaloids" = "#ec6669", "Phenolic acids" = "#f39668", "Fatty acids" = "#bd647d", 
              "Carbonyl compounds" = "#832c9b", "Monosaccharides" = "#5f70a8", "Cholines" = "#d4a34b", "TCA acids" = "#4da457", 
              "Azoles" = "#d16fa8", "Organonitrogen compounds" = "#78b646", "Carboxylic acids" = "#3b738f", 
              "Short-chain acids" = "#a36f2d", "Sulfonic acids" = "#b54f82", "Amines" = "#6e2b57", "Amine oxides" = "#c97f2e", 
              "Organic carbonic acids" = "#8e4555", "Hydrocarbons" = "#4f6d32", "Phenylpropanoids" = "#395262", 
              "Primary alcohols" = "#9c572c", "Alcohols and polyols" = "#5a9d79", "Phenols" = "#d09636", 
              "Fatty amines" = "#7d3b65", "Pyrimidines" = "#1b69a1", "Purines" = "#48a0af")
  
  

# Metabolites that are unknown to RefMet - try more standardized name
t_results$variable[t_results$variable=="Acetate.mM."] <- "Acetate"
t_results$variable[t_results$variable=="Methyl.2.oxovalerate"] <- "3-Methyl-2-oxovaleric acid"
t_results$variable[t_results$variable=="Aminopentanoate"] <- "5-Aminopentanoate"
t_results$variable[t_results$variable=="TMA..N.oxide"] <- "Trimethylamine N-oxide"
t_results$variable[t_results$variable=="X2.Hydroxyisobutyrate"] <- "2-hydroxybutyric acid"
t_results$variable[t_results$variable=="X3.Hydroxyisobutyrate"] <- "3-hydroxybutyric acid"

t_results$variable <- gsub("^X\\d+\\.", "", t_results$variable)

# Replace all "-" entries with "unknown" in merged_data dataframe
RefMet_mapped <- RefMet_mapped %>%
  mutate_all(~replace(., . == "-", "unknown"))


merged_data <- merge(RefMet_mapped, t_results, by.x = "Input.name", by.y = "variable", all.x = F, all.y = TRUE)



# Create a new dataframe with the same column names
new_df <- data.frame(Sub.class = merged_data[,class],
                     p_M2 = rep(NA, nrow(merged_data)),
                     p_M4 = rep(NA, nrow(merged_data)),
                     p_M6 = rep(NA, nrow(merged_data)))

colnames(new_df)[1] <- "Overall"


# Fill in 'yes' where conditions are met
new_df$p_M2[merged_data$p_M2 < 0.05] <- new_df$Overall[merged_data$p_M2 < 0.05]
new_df$p_M4[merged_data$p_M4 < 0.05] <- new_df$Overall[merged_data$p_M4 < 0.05]
new_df$p_M6[merged_data$p_M6 < 0.05] <- new_df$Overall[merged_data$p_M6 < 0.05]

# Do not plot reference
if (is.null(plot.reference)) {
  new_df[,1]<-NULL}


tmp<- new_df %>%
  tidyr::pivot_longer(cols = everything(), names_to = "time", values_to = "subname")

tmp <- na.omit(tmp)

tmp <- as.data.frame(table(tmp$time, tmp$subname))

tmp <- tmp %>%
  group_by(Var1) %>%
  mutate(sum = sum(Freq)) %>%
  ungroup()

tmp <- tmp |> 
  mutate(prop = Freq / sum,
         n = sum) %>%
  group_by(Var1) |> 
  arrange(desc(Var2)) |> 
  mutate(ypos = cumsum(prop) - 0.5*prop,
         test = cumsum(prop) ) |> 
  ungroup()

# New facet label names 
new <- c(paste0("<b>Total</b><br>(n=",tmp$n[1],")"), 
         paste0("<b>M2</b><br>_p<0.05_<br>(n=",tmp$n[2],")"),
         paste0("<b>M4</b><br>_p<0.05_<br>(n=",tmp$n[3],")"),
         paste0("<b>M6</b><br>_p<0.05_<br>(n=",tmp$n[4],")"))
names(new) <- c("Overall", "p_M2","p_M4","p_M6")


p <- tmp |> 
  dplyr::filter(Freq != 0) |> 
  dplyr::group_by(Var1) |> 
  dplyr::mutate(label = if_else(Freq/n >= 0.3, paste0("<span style='color:white'>", Freq, "/", n, "<br></span>"), NA)) |> 
  dplyr::ungroup() |> 
  ggplot(aes(x = Var1,
             y= prop,
             fill = Var2)) +
  coord_flip() +
  geom_bar(stat = 'identity',
           width = 0.8) +
  facet_wrap(Var1~.,
             scales = 'free_y', ncol=1, strip.position="left", 
             labeller = labeller(Var1 = new)) +
  geom_richtext(aes(y = ypos,
                    label = label),
                size = 2.3,
                label.color = NA,
                fill = NA,
                show.legend = F) +
  theme_void() +
  theme(legend.position = "right",
        strip.text=element_markdown()) +
  scale_colour_manual(values = cols,
                      aesthetics = "fill",
                      name = legend.name) +
  theme(legend.key.height= unit(0.5, 'cm'),
        legend.key.width= unit(0.5, 'cm'),
        legend.text = element_text(size=8))


return(p)

}