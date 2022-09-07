## Stats for FACS

library(dunn.test)
library(dplyr)

## Formating and combining dfs - need to have one df with columns in desired order (can be an empty df just specifying cols)

# Format function
## metadata needs to have rownames same as df


organiser <- function(df, desired.rows, desired.cols, df.order, metadata, v1=NULL, 
                      v2=NULL, v3=NULL){
  
 
  sample.names <- data.frame(SampleID = rownames(df)[desired.rows])
  df <- df[desired.rows, desired.cols]
  rownames(df) <- sample.names
  metadata <- metadata[match(rownames(df), rownames(metadata))]
  
  if(!(is.null(v1))){df[,v1] <- metadata[,v1]}
  if(!(is.null(v2))){df[,v2] <- metadata[,v2]}   
  if(!(is.null(v3))){df[,v3] <- metadata[,v3]}
  
  df <- setcolorder(df, names(df.order))
  
  return(df)
  }

## Running the formatting on a list of dfs
data.organiser <- function(data.frame.list, desired.rows, desired.cols, df.order, metadata, v1=NULL, v2=NULL, v3=NULL){
  
  temp.df.list <- list()
for(i in 1:length(data.frame.list)){
  spruced.df <- organiser(data.frame.list[[i]], desired.rows, desired.cols, df.order, metadata, v1, v2, v3)
  temp.df.list[[i]] <- spruced.df}
  fixed.df <- bind_rows(temp.df.list)
  return(fixed.df)
}

## Convert frequencies to cell counts

freq2counts <- function(freq.df, desired.cell.cols, cell.counts){
  cell.counts <-  cell.counts[match(rownames(freq.df), rownames(cell.counts)),] 
  freq.df <- freq.df[,desired.cell.cols]
  for(i in 1:ncol(freq.df)){
    freq.df[,i] <- freq.df[,i] * cell.counts
  }
  
  return(freq.df)
}

## Kruskal wallis with dunns per column in a FACS data frame

## Option for subsetting data e.g. based on tissue

## Requires metadata with iv - 
## For FACS I used paste0 of Condition + day for individual samples obtained on different days

## The actual dunn test

library(dunn.test)

dunn <- function(col, metadata, independent.variable, adj){
  iv <- unlist(metadata[, independent.variable])
  kw.result <- kruskal.test(unlist(col) ~ as.factor(iv))
  posthoc <- dunn.test(unlist(col), as.factor(iv), method = adj, list = TRUE)
  posthoc.res <- data.frame(t(data.frame(posthoc$P.adjusted)))
  
  names(posthoc.res) <- posthoc$comparisons
  
  kw <-  data.frame(kw.p = kw.result$p.value)
  final.res <- cbind(kw, posthoc.res)
  return(final.res)
}

### Running the dunn test per column in df

run.dunn <- function(df, subset=NULL, subsetcol, metadata, independent.variable, testnum, adjust){
  if (!(is.null(subset))){
    
    metadata <- metadata[metadata[, subsetcol] == subset, ]}
  df <- df[rownames(df) %in% rownames(metadata),]
  df <- df[,1:testnum]
  aov.results=list()
  for(i in 1:ncol(df)){
    result <- dunn(df[,i], metadata, independent.variable, adjust)
    aov.results[[i]] <- as.data.frame(result)
    
  }
  aov.results <- bind_rows(aov.results)
  rownames(aov.results) <- colnames(df)
  return(aov.results)
  
}

## starry.eyed
################
# Coverts P values given above into stars

### Gives back a new df with stars instead of raw p values - have this as new df rather than overwriting original stats
## Answer given by Rui Barradas 
## https://stackoverflow.com/questions/57231854/substitute-p-values-for-stars-in-data-frame-in-r

stars.pval <- function(x){
  stars <- c("****", "***", "**", "*", "ns")
  vec <- c(0, 0.0001, 0.001, 0.01, 0.05, 1)
  i <- findInterval(x, vec)
  stars[i] }

## Function for iterating over df by column - use to turn your stats results table into stars :) 

starry.eyed <- function(df){
  tempdf <- list()
  for(i in 1:ncol(df)){
    stars <- stars.pval(df[,i])
    tempdf[[i]] <- as.data.frame(stars)
  }
  
  stars.res <- bind_cols(tempdf)
  colnames(stars.res) <- colnames(df)
  rownames(stars.res) <- rownames(df)
  return(stars.res)
  
}

########## Wilcox.test 
########## E.g. comparison of 1 cell type between 2 conditions only

wil.test <- function(col, metadata, independent.variable){
  iv <- unlist(metadata[, independent.variable])
  wil.result <- wilcox.test(unlist(col) ~ as.factor(iv))
  result <- data.frame(wilcox.p = wil.result$p.value)
  return(result)}


########### Run over a data frame by column

run.wil.test <- function(df, subset=NULL, subsetcol, metadata, independent.variable, col4test){
  if (!(is.null(subset))){ metadata <- metadata[metadata[, subsetcol] == subset, ]}
  df <- df[rownames(df) %in% rownames(metadata),]
  df <- df[,col4test]
  wil.results = list()
  for(i in 1:ncol(df)){
    result <- wil.test(df[,i], metadata, independent.variable)
    wil.results[[i]] <- as.data.frame(result)
    
  }
  wil.results <- bind_rows(wil.results)
  wil.results$padj <- p.adjust(wil.results$wilcox.p, method = "BH")
  wil.results$padj.star <- stars.pval(wil.results$padj)
  if (!(is.null(subset))){wil.results[,subsetcol] <- rep(subset, nrow(wil.results))}
  rownames(wil.results) <- colnames(df)
  return(wil.results)}











