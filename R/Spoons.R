###############################################
###############################################
################## Spoons #####################
###############################################
###############################################
# For scooping
###############################################
# Includes:
# TOIPlucker, BoxPlotTOI, Tungstify, FACSPlotter, FACSPlotterXtreme
###############################################

library(tidyr)
library(ggplot2)
library(data.table)

###############################################
###############################################
# TOIPlucker
###############################################
# Pull out "thing of interest" from a normalised
# matrix e.g. gene name from transcript abundance data
# 
# Note: rownames of metadata have to be same as colnames of normalised.matrix
# Note 2: "thing of interest" has to be in rownames of normalised.matrix
################################################

TOIPlucker <- function(TOI, normalised.matrix, metadata, group = "Condition"){
TOI_Abundance <- data.frame(normalised.matrix[TOI,])
TOI_Abundance <- TOI_Abundance[,rownames(metadata)]
TOI_Abundance$TOI <- rownames(TOI_Abundance)
TOI_L <- gather(TOI_Abundance, "Samples", "Value", -TOI)
TOI_L$group <- metadata[TOI_L$Samples, group]

return(TOI_L)

}

###############################################
###############################################
# BoxPlotTOI
###############################################
# Standard boxplot of your "thing of interest"
# y axis will be normalised value 
# x will be group you specify
###############################################

BoxPlotTOI <- function(TOI_L){

Plot <- ggplot(TOI_L, aes(x = group, y = Value, color = TOI)) +
geom_boxplot(outlier.alpha = 0) +
geom_point(position = position_jitterdodge(jitter.height = 0, jitter.width = 0.1), alpha = 0.5) +
theme_classic() +
  theme(text = element_text(size=20), axis.text.x = element_text(size = 20, color = "black", 
                                                                 margin = margin(t=0, r=0, l=0, b=5)), axis.text.y = element_text(size = 20, color = "black",  margin = margin(t=0, r=0, l=5, b=0))) 
  

return(Plot)
}

##############################################
##############################################
# Tungstify
##############################################
# Make tungstate graphs the right colours with big text size for axis
#
#For graph consistency at FP meetings
#
# Assumes tungstate condition is specified by color = Condition -
# in an order where control comes first
#
# Probs abit redundant now as have the FACSPlotter functions 
# but like the name too much to delete
##############################################

Tungstify <- function(Tungstate.Graph){
  
Plot <-  Tungstate.Graph +
  scale_color_manual(values = c("#6d6d6dff", "#7900f5ff")) +
  theme_classic() +
  theme(text = element_text(size=20), axis.text.x = element_text(size = 20, color = "black", 
                                                                 margin = margin(t=0, r=0, l=0, b=5)), axis.text.y = element_text(size = 20, color = "black",  margin = margin(t=0, r=0, l=5, b=0))) 
  
return(Plot)  
}

############################################
############################################
# FACSPlotter
############################################
# FACSPlotter: Vanilla
############################################
# The OG function to plot a singular boxplot showing one aspect of your FACS Data
############################################

FACSPlotter <- function(FACS.data, x.axis, y.axis, colour.by, colour.by.title = "Condition", x.title = "Day", y.title = "FOXP3 T Cells (%)", legend.position = "right", facet = FALSE, facet.by = "none"){
  Plot <- ggplot(FACS.data, aes_string(x = x.axis, y = y.axis, color = colour.by)) +
    geom_boxplot(outlier.alpha = 0) +
    geom_point(position = position_jitterdodge(jitter.height = 0, jitter.width = 0.1), alpha = 0.5) +
    scale_color_manual(values = c("#6d6d6dff", "#7900f5ff"), name = colour.by.title) +
    theme_classic() +
    theme(legend.position = legend.position, text = element_text(size=20), axis.text.x = element_text(size = 20, color = "black", 
                                                                   margin = margin(t=0, r=0, l=0, b=5)), axis.text.y = element_text(size = 20, color = "black",  margin = margin(t=0, r=0, l=5, b=0))) +
    labs(x = x.title, y = y.title) + 
    expand_limits(y = 0)  #Makes sure each graph starts from 0
  
  if(facet == TRUE){Plot <- Plot + facet_wrap(facet.by)}
 
  return(Plot)
  
}

#############################################
#############################################
# FACSPlotterXtreme
#############################################
# FACSPlotter: Extreme Edition
#############################################
# Churns out multiple FACS Plots for variables that are columns in a FACS data dataframe 
# e.g. multiple cells of interest
#
# Have the option to save the plots if you would like
#############################################

FACSPlotterXtreme <- function(FACS.data, x.axis, y.axis.as.list, colour.by, 
                              colour.by.title = "Condition", x.title = "Day", y.title.as.list = y.titles.as.list,
                              legend.position = "right", 
                              file.name = "Spleen", file.width = 18.5, file.height = 16, save = FALSE, facet = FALSE, facet.by = "none"){
  
  plot_list = list()
  var_list = list(cell_plot = y.axis.as.list, 
                  y_label = y.title.as.list)
  for (i in 1:length(var_list$cell_plot)){
    
    Plot <- FACSPlotter(FACS.data, x.axis, var_list$cell_plot[i], colour.by, colour.by.title = colour.by.title, y.title = var_list$y_label[i], x.title = x.title, legend.position = legend.position, facet = facet, facet.by = facet.by)
      plot_list[[i]] = Plot
      if(save == TRUE){ggsave(Plot, file=paste0(file.name, var_list$cell_plot[i], ".png"), width = file.width, height = file.height, units = "cm")}
  }
  
  return(plot_list)
}

##############################################
##############################################
# TOIbySpeciesPlotter

TOIbySpecPlotter <- function(TOI_L, x.axis, y.axis, fill.by, facet = FALSE, facet.by = "none", save = FALSE, file.name = "none", file.width = 18.5, file.height = 16, x.title = "X", y.title = "Y"){
  
  Plot <- ggplot(TOI_L, aes_string(x = x.axis, y = y.axis, fill = fill.by)) + 
    geom_col() +
    theme_classic() +
    scale_fill_manual(values = c("steelblue2", "darkseagreen", "palevioletred2", "orange1", "plum3", "lightblue3", "sandybrown", "lightpink", "darkolivegreen2", "lightsteelblue1", "lightgrey")) +
   labs(x = x.title, y = y.title)
  
  
  
  if(facet == TRUE){Plot <- Plot + facet_wrap(facet.by)}
  if(save == TRUE){ggsave(Plot, file=paste0(file.name,".png"), width = file.width, height = file.height, units = "cm")}
  return(Plot)
}
  

################################
# TOIbySpecUltimate

TOIbySpecUlt <- function(infiles, path.to.TOI.file, infile.dir, metadata, total.row.means = TRUE, total.row.means.alternative = 1:24,
                          group.1.name = "Condition", group.2 = TRUE, group.2.name = "Day", group.3 = TRUE, group.3.name = "Mouse", x.axis.to.plot = "Mouse", 
                         y.axis.to.plot = "Rel_Abun", fill.to.plot = "Bacteria", facet = TRUE, facet.by = "Day", y.axis.label = "Relative Abundance (%)", x.axis.label = "Mouse"){
  
  
  GrobsList <- list()
  
  
  for(i in 1:length(infiles)){
    
    Gene.name <- gsub(".tsv", "", unlist(strsplit(infiles[i], "_"))[2])
    Y.name <- paste0(Gene.name, y.axis.label)
    path.to.infiles <- paste0(infile.dir, infiles[i])
    
    Spec.df <- TOIbySpecies(path.to.TOI.file, path.to.infiles, metadata, total.row.means = total.row.means, total.row.means.alternative = total.row.means.alternative, 
                            group.1.name = group.1.name, group.2 = group.2, group.2.name = group.2.name, group.3 = group.3, group.3.name = group.3.name)  
    Spec.plot <- TOIbySpecPlotter(Spec.df, x.axis.to.plot, y.axis.to.plot, fill.to.plot, facet = facet, facet.by = facet.by, y.title = Y.name, x.title = x.axis.label)
    GrobsList[[i]] <- Spec.plot 
    
    }
  
  return(GrobsList)
  }




