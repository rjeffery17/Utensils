library(data.table)
library(DESeq2)

#######################################
#######################################
#######################################
#Calculate abundance from count data
#NB/ gene name or w/e cannot be a column in count.matrix
#Import files with row.names = 1

CalcAbun <- function(count.matrix){

Abun_Matrix <- as.data.frame(t((t(count.matrix)/colSums(count.matrix))*100))
  
return(Abun_Matrix)
}

######################################
######################################
# Add readable descriptions to things of interest e.g. multiple KEGGs grouped by function
# Need $Description in TOI data - filled with NAs
# TOI.group = list of TOI names which can be grouped under one description

AddDescript <- function(TOIs.data, TOI.group, description="tup"){

  TOIs.data$Description <- ifelse(TOIs.data$TOI %in% TOI.group, description, TOIs.data$Description)
  return(TOIs.data)
}


#######################################
#######################################
#Import, remove mean and SD and % signs from raw flowjo export files

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
#TOIbySpecies
########################################
# Create table showing TOI abundance over time with species contribution to TOI

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
## Multi-Omic Longituinal DESEq
###############################
# Run DESeq on TOI data from two different time points on two different omic datasets


LongDESeq <- function(DNA.metadata, RNA.metadata, day.1, day.2, DNA.count.data, RNA.count.data){
  DNA_metadata <- DNA.metadata[DNA.metadata$Day %in% c(day.1, day.2),]
  RNA_metadata <- RNA.metadata[RNA.metadata$Day %in% c(day.1, day.2),]
  DNA_count_data <- DNA.count.data[, rownames(DNA_metadata)]
  RNA_count_data <- RNA.count.data[,rownames(RNA_metadata)]
  
  # Do the DDS for DNA
  DNA_metadata$Day <- factor(DNA_metadata$Day)
  DNA_metadata$Mouse <- factor(DNA_metadata$Mouse)
  DNA_count_data_dds <- DESeqDataSetFromMatrix(countData = DNA_count_data,
                                           colData = DNA_metadata, 
                                           design = ~ Mouse + Day)
  
  DNA_count_data_dds$Day <- relevel(DNA_count_data_dds$Day, ref = "-3")
  
  DNA_count_data_dds <- DESeq(DNA_count_data_dds, test = "LRT", fitType = "local", reduced = ~ Mouse)
  DNA_dds_results <- results(DNA_count_data_dds) 
 
  #Same for RNA
  RNA_metadata$Day <- factor(RNA_metadata$Day)
  RNA_metadata$Mouse <- factor(RNA_metadata$Mouse)
  RNA_count_data_dds <- DESeqDataSetFromMatrix(countData = RNA_count_data,
                                               colData = RNA_metadata, 
                                               design = ~ Mouse + Day)
  
  RNA_count_data_dds$Day <- relevel(RNA_count_data_dds$Day, ref = "-3")
  
  RNA_count_data_dds <- DESeq(RNA_count_data_dds, test = "LRT", fitType = "local", reduced = ~ Mouse)
  RNA_dds_results <- results(RNA_count_data_dds) 
  
  RNA_DNAK <- data.frame(Gene=rownames(DNA_count_data), RNA_FC=RNA_dds_results$log2FoldChange, DNA_FC=DNA_dds_results$log2FoldChange, RNA_padj=RNA_dds_results$padj, DNA_padj=DNA_dds_results$padj)
  
  return(RNA_DNAK)
  
}

################################
################################
#### Gene Plucker
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



