library(data.table)
library(DESeq2)

#######################################
#######################################
# CalcAbun
#######################################
# Calculate abundance from count data
# NB/ gene name or w/e cannot be a column in count.matrix
# Import files with row.names = 1

CalcAbun <- function(count.matrix){

Abun_Matrix <- as.data.frame(t((t(count.matrix)/colSums(count.matrix))*100))
  
return(Abun_Matrix)
}

######################################
######################################
# AddDescript
######################################
# Add readable descriptions to things of interest e.g. multiple KEGGs grouped by function
# Need $Description in TOI data - filled with NAs
# TOI.group = list of TOI names which can be grouped under one description
######################################
# Redundant - use DESCRIPTOR2000
######################################

AddDescript <- function(TOIs.data, TOI.group, description="tup"){

  TOIs.data$Description <- ifelse(TOIs.data$TOI %in% TOI.group, description, TOIs.data$Description)
  return(TOIs.data)
}

#######################################
#######################################
# FACSSprucer
#######################################
#Import, remove mean and SD and % signs from raw flowjo export files
#######################################

FACSSprucer <- function(path.to.file, sample.number = 10, conditions = Conditions, days = Days, tissue = Tissue){
  Data <- read.csv(path.to.file, stringsAsFactors = FALSE, row.names = 1)
  Data <-  Data[1:sample.number,]
  Data$Sample <- row.names(Data)
  Data$Condition <- conditions
  Data$Day <- days
  Data$Tissue <- tissue
 
  PercentEraser <- function(vector.of.strings){
    
    newvector <- gsub(" %","", vector.of.strings)
    return(newvector)
  }
  
  Data[,!colnames(Data) %in% c("Sample", "Condition", "Tissue", "Day")] <- as.data.frame(apply(Data[,!colnames(Data) %in% c("Sample", "Condition", "Tissue", "Day")], 2, PercentEraser))
  Data[,!colnames(Data) %in% c("Sample", "Condition", "Tissue", "Day")] <- apply(Data[,!colnames(Data) %in% c("Sample", "Condition", "Tissue", "Day")], 2, function(x) as.numeric(as.character(x)))
  return(data.frame(Data))
}

########################################
########################################
# TOIbySpecies
########################################
# Create table showing TOI abundance over time with species contribution to TOI
#########################################

TOIbySpecies <- function(path.to.TOI.file, path.to.TOI.by.species.file, metadata, total.row.means = TRUE, total.row.means.alternative = 1:24,
                         group.1.name = "Condition", group.2 = TRUE, group.2.name = "Day", group.3 = TRUE, group.3.name = "Mouse"){
  TOI <- read.csv(path.to.TOI.by.species.file, sep = "\t", row.names = 1, stringsAsFactors = FALSE)
  # Subset TOI according to samples in metadata 
  TOI <- TOI[ , rownames(metadata)]
  TOI_Rel_Abun <- read.csv(path.to.TOI.file, sep = "\t", row.names = 1, stringsAsFactors = FALSE)
  TOI_Rel_Abun <- TOI_Rel_Abun[, rownames(metadata)]
  
  #Order Files according to metadata
  TOI <- setcolorder(TOI, rownames(metadata))
  TOI_Rel_Abun <- setcolorder(TOI_Rel_Abun, rownames(metadata))
  
  #Get each KEGG by Species relative to the abundance of all counts
  TOI_Con <- as.data.frame(t((t(TOI)/colSums(TOI_Rel_Abun))*100))
  
  #Order based on rowmeans or select sample row means
  TOI_Con <- TOI_Con[order(rowMeans(TOI_Con), decreasing = TRUE),]
  if(total.row.means == FALSE){TOI_Con <- TOI_Con[order(rowMeans(TOI_Con[,total.row.means.alternative]), decreasing = TRUE),]}
  TOI_Con_Top10 <- TOI_Con[1:10,]
  TOI_Con_Top10 <- rbind(TOI_Con_Top10, colSums(TOI_Con[11:nrow(TOI_Con),]))
  rownames(TOI_Con_Top10)[11] <- "Other"
  TOI_Con_Top10$Bacteria <- unlist(strsplit(rownames(TOI_Con_Top10), "f__"))[c(2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 21)]
  
  TOI_L <- gather(TOI_Con_Top10, 'Samples', 'Rel_Abun', -Bacteria)
  TOI_L$group1 <- metadata[TOI_L$Samples , group.1.name]
  if(group.2 == TRUE){TOI_L$group2 <- metadata[TOI_L$Samples, group.2.name]}
  if(group.3 == TRUE){TOI_L$group3 <- metadata[TOI_L$Samples, group.3.name]}
  TOI_L <- setnames(TOI_L, "group1", group.1.name)
  if(group.2 == TRUE){TOI_L <- setnames(TOI_L, "group2", group.2.name)}
  if(group.3 == TRUE){TOI_L <- setnames(TOI_L, "group3", group.3.name)}
  
 #took this out becuase if you have 2 unassigned in the top 10 it throws a wobbler :(  
  #TOI_L$Bacteria <- factor(TOI_L$Bacteria, levels = TOI_Con_Top10$Bacteria)
  return(TOI_L)
}

###############################
###############################
## Multi-Omic Longituinal DESEq
###############################
# Run DESeq on TOI data from two different time points on two different omic datasets
###############################


LongDESeq <- function(RNA_MetaData, RNA_KEGG, DNA_MetaData, DNA_KEGG, day.1, day.2){
  
  
  RNA_MetaData <- RNA_MetaData[RNA_MetaData$Day %in% c(day.1, day.2),]
  DNA_MetaData <- DNA_MetaData[DNA_MetaData$Day %in% c(day.1, day.2),]
  RNA_KEGG <- RNA_KEGG[, rownames(RNA_MetaData)]
  DNA_KEGG <- DNA_KEGG[, rownames(DNA_MetaData)]
  
  
  RNA_MetaData$Day <- as.factor(RNA_MetaData$Day)
  DNA_MetaData$Day <- as.factor(DNA_MetaData$Day)
  
  RNAKEGG.dds <- DESeqDataSetFromMatrix(countData = RNA_KEGG,
                                        colData = RNA_MetaData, 
                                        design = ~ Mouse + Day)
  
  DNAKEGG.dds <- DESeqDataSetFromMatrix(countData = DNA_KEGG,
                                        colData = DNA_MetaData, 
                                        design = ~ Mouse + Day)
  
  RNAKEGG.dds$Day <- relevel(RNAKEGG.dds$Day, ref = "-3")
  DNAKEGG.dds$Day <- relevel(DNAKEGG.dds$Day, ref = "-3")
  
  RNAKEGG.dds <- DESeq(RNAKEGG.dds, test = "LRT", fitType = "local", reduced = ~Mouse) 
  DNAKEGG.dds <- DESeq(DNAKEGG.dds, test = "LRT", fitType = "local", reduced = ~Mouse) 
  
  RNAKEGG.res <- results(RNAKEGG.dds)
  DNAKEGG.res <- results(DNAKEGG.dds)
  RNA_DNAK <- data.frame(Gene=rownames(RNA_KEGG), RNA_FC=RNAKEGG.res$log2FoldChange, DNA_FC=DNAKEGG.res$log2FoldChange, RNA_padj=RNAKEGG.res$padj, DNA_padj=DNAKEGG.res$padj)
  return(RNA_DNAK)
}

################################
################################
#### Gene Plucker
################################
# Pluck out multiple genes with a common signature e.g. nap will return napA, napB etc.
################################

GenePlucker <- function(Gene_List, DF, variable = "Gene"){
  #Create an empty vector to store results
  res <- c()
  for (x in Gene_List){
    genes <- DF[grep(x, DF[, variable]),]
    genes <- as.character(genes[, variable])
    res <- append(res, genes)
  }
  
  return(res)
}

################################
################################
#  TOIPluckerXtreme
################################

TOIPluckerXtreme <- function(TOI.list, df, TOI.variable.name.in.df){
  Data <- df[df[,TOI.variable.name.in.df] %in% TOI.list, ]
  return(Data)
}

################################
################################
# DESCRIPTOR2000
################################
# Adds descriptions based on a list where the key is the readable decription 
# name and the values are the things to assign to that description
################################

DESCRIPTOR2000 <- function(descriptions, data.table){
  data.table$Description <- rep(NA, nrow(data.table))
  for(i in 1:length(descriptions)){
    data.table$Description <- ifelse(data.table$Gene %in% descriptions[[i]], names(descriptions)[i], data.table$Description)
  } 
  return(data.table)
}

################################
################################
# fancy_scientific
################################
# Scientific labelling of log axis
###############################

fancy_scientific <- function(l) {
  # turn in to character string in scientific notation
  l <- format(l, scientific = TRUE)
  # quote the part before the exponent to keep all the digits
  l <- gsub("^(.*)e", "'\\1'e", l)
  # turn the 'e+' into plotmath format
  l <- gsub("e", "%*%10^", l)
  # return this as an expression
  parse(text=l)
}


