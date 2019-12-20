
library(tidyr)

###############################################
###############################################


CreateTOI <- function(TOI, normalised.matrix, metadata, group = "Condition"){
TOI_Abundance <- data.frame(normalised.matrix[TOI,])
TOI_Abundance <- TOI_Abundance[,rownames(metadata)]
TOI_Abundance$TOI <- rownames(TOI_Abundance)
TOI_L <- gather(TOI_Abundance, "Samples", "Value", -TOI)
TOI_L$group <- metadata[TOI_L$Samples, group]

return(TOI_L)

}


###############################################
###############################################

BoxPlotTOI <- function(TOI_L){

Plot <- ggplot(TOI_L, aes(x = group, y = Value)) +
geom_boxplot(outlier.alpha = 0) +
geom_jitter(alpha = 0.5) +
theme_classic() +
  ggtitle(TOI)

return(Plot)
}