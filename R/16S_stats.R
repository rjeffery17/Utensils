#### Data Organisation and stats for 16S 

library(tidyr)
library(reshape2)


######################################
### Order based on rowsums
######################################

order.rows <- function(abundance.matrix){
  abundance.matrix <- abundance.matrix[order(rowSums(abundance.matrix), decreasing = TRUE),]
}

#######################################
#######################################
# CalcAbun
#######################################
# Calculate abundance from count data
# NB/ gene name or w/e cannot be a column in count.matrix
# Import files with row.names = 1

CalcAbun <- function(count.matrix){
  
  Abun_Matrix <- as.data.frame(t((t(count.matrix)/colSums(count.matrix))*100))
  Abun_Matrix <- order.rows(Abun_Matrix)
  return(Abun_Matrix)
}

#######################################
#######################################


## Takes bacteria from rownames of dada2 output and returns the genus name (or wherever you choose to split)
## Also removes incertae sedis 

genus.name <- function(df, split.point = "g__"){
  
  bug.names <- unlist(strsplit(row.names(df), split.point))
  bug.names <- bug.names[seq(2, length(bug.names), 2)]
  bug.names <- gsub(" incertae sedis", "", bug.names)
  return(bug.names)
  
}

## Make data long - adds bacteria as output of genus.name function

make.long <- function(df, meta, iv, time.point, var3=NULL, rename="g__"){
  if(!is.null(rename)){df$Bacteria <- genus.name(df, rename)}
  if(is.null(rename)){df$Bacteria <- rownames(df)}
  df <- gather(df, "Sample", "Rel_Abun", -Bacteria)
  meta <- meta[match(df$Sample, rownames(meta)),]
  df[,iv] <- meta[,iv]
  df[,time.point] <- meta[,time.point]
  df[,time.point] <- as.factor(df[,time.point])
  if(!is.null(var3)){df[, var3] <- meta[,var3]}
  return(df) }


################ Add strain names to genus #############
########################################################

bacterial.strains <- function(df, strain.dict){
  df$bacteria.strain <- rep(NA, nrow(df))
  for(i in 1:length(strain.dict)){
    df$bacteria.strain <- ifelse(df$Bacteria %in% names(strain.dict[i]), 
                                paste0(df$Bacteria, " (", strain.dict[[i]], ")"), df$bacteria.strain)
  }
  
  df$bacteria.strain <- ifelse(is.na(df$bacteria.strain), as.character(df$Bacteria), df$bacteria.strain)
  return(df)
}

###############################################################################
################### Turn longitudinal readings into average readings per sample
################# E.g. Used by me to get average of naive stool readings to then compare
################## stool with caecum, colon, SI
#############################################################################
####### Bits hard-coded = SampleID has to be column name with sample identifier in
######################### bacteria.strain has to be column name with bacteria in it

#### Function to get average abundance of each bacteria per sample

avg.rel.abun <- function(bug.list, df){
  df.return = list()
  for(i in 1:length(bug.list)){
    tmp.bug <- df[df$bacteria.strain == bug.list[i],]
    tmp.bug.out <- data.frame(tmp.bug[1, c(1:2, 4:8)])
    tmp.bug.out$Avg_Rel_Abun <- mean(tmp.bug$Rel_Abun)
    df.return[[i]] <- tmp.bug.out
  }
  
  avg <- bind_rows(df.return)
  return(avg)  
}

############# Function to iterate through samples in a data frame

run.per.sample.avg.abun <- function(SampleID.list=sample.IDs, bug.list=mm12.list, df){
  avgs <- list()
  for(i in 1:length(SampleID.list)){
    single.sample <- df[df$SampleID == SampleID.list[i],]
    sample.avg <- avg.rel.abun(bug.list, single.sample)
    avgs[[i]] <- sample.avg
  }
  
  avgs <- bind_rows(avgs)
  return(avgs)
  
}


########################################################
########################################################
######### lmer Longitudinal Analysis ###################
########################################################
########################################################

library(lmerTest)
library(emmeans)
library(dplyr)

#### Convert p values to stars - also in FACS_stats
## Answer given by Rui Barradas 
## https://stackoverflow.com/questions/57231854/substitute-p-values-for-stars-in-data-frame-in-r

stars.pval <- function(x){
  stars <- c("****", "***", "**", "*", "ns")
  vec <- c(0, 0.0001, 0.001, 0.01, 0.05, 1)
  i <- findInterval(x, vec)
  stars[i] }

## Function for iterating over df by column - use to turn your dunn.result table into stars :) 

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

########################################################
### The actual lmer testing 
### NB metadata must have rownames matching counts colnames

lmr.emeans <- function(row, metadata, independent.variable="Condition", sample.variable="SampleID", time.variable="Day", 
                       emeans.factor, emeans.by.factor){
  
  bug.name <- rownames(row)
  bug.df <- data.frame(Sample=names(row), Rel_Abun=unlist(row))
  
  # bug.df$Bacteria <- bug.name
  
  # bug.df <- gather(bug.df, "Sample", "Rel_Abun", -Bacteria)
  
  metadata <- metadata[match(bug.df$Sample, rownames(metadata)), ]
  
  bug.df$time <- metadata[,time.variable]
  bug.df$independent.variable <- metadata[,independent.variable]
  bug.df$sample <- metadata[,sample.variable]
  
  bug.df$time <- as.factor(bug.df$time)
  bug.df$independent.variable <- as.factor(bug.df$independent.variable)
  
  
  ## lmr Model
  lmr.res <- lmer(Rel_Abun ~ independent.variable * time + (1 | sample), data = bug.df)
  
  # Anova on lmr
  summary <- anova(lmr.res)
  anova.res <- t(data.frame(summary$P))
  rownames(anova.res) <- bug.name
  colnames(anova.res) <- c(independent.variable, time.variable, paste0(independent.variable, "*" ,time.variable))
  anova.res <- as.data.frame(anova.res)
  
  # Pairwise on lmr
  emmeans.res <- emmeans(lmr.res, emeans.factor, by = emeans.by.factor)
  results <- pairs(emmeans.res)
  results <- data.frame(results)

  
  bacteria <- data.frame(Bacteria = rep(bug.name, nrow(results)))
  test.complete <- cbind(bacteria, results)
  
  overall.results <- list(anova.res, test.complete)
  names(overall.results) <- c("lmr.output", "pairwise.emmeans.output")
  
  return(overall.results)
  
}

##############################################################
#### Execute on multiple rows in a counts table ##############

run.lmer.emeans <- function(normalised.counts, metadata, independent.variable="Condition", sample.variable="SampleID", time.variable="Day", emeans.factor="independent.variable", emeans.by.factor="time"){
  
  anova.results <- list()
  pairwise.results <- list()
  
  for (i in 1:nrow(normalised.counts)){
    result <- lmr.emeans(normalised.counts[i,], metadata, independent.variable, sample.variable, time.variable, emeans.factor, emeans.by.factor)
    anova.results[[i]] <- as.data.frame(result$lmr.output)
    pairwise.results[[i]] <- as.data.frame(result$pairwise.emmeans.output)
  }
  
  anova.results <- bind_rows(anova.results)
  rownames(anova.results) <- rownames(normalised.counts)
  anova.results[,4] <- p.adjust(anova.results[,1], method = "BH")
  anova.results[,5] <- p.adjust(anova.results[,2], method = "BH")
  anova.results[,6] <- p.adjust(anova.results[,3], method = "BH")
  colnames(anova.results)[4:6] <- c(paste0(independent.variable, ".padj"), paste0(time.variable, ".padj"), paste0(independent.variable, "*" ,time.variable, ".padj"))
  
  
  pairwise.results <- bind_rows(pairwise.results)
  
  ## P value adjusted in Tukey method - not sure if BH required on top of this
  pairwise.results$padj <- p.adjust(pairwise.results$p.value, method = "BH")
  pairwise.results$p.star <- stars.pval(pairwise.results$padj)
  
  
  all.res <- list(anova.results, pairwise.results)
  names(all.res) <- c("lmr.aov.out", "emmeans.out")
  return(all.res) }

###############################################################
###############################################################
## Kruskal Wallis with Dunns for End-Point Data
###############################################################

library(dunn.test)

dunn <- function(row, metadata, independent.variable, adj){
  
  asv <- rownames(row)
  iv <- unlist(metadata[,independent.variable])
  
  kw.result <- kruskal.test(unlist(row) ~ as.factor(iv))
  posthoc <- dunn.test(unlist(row), as.factor(iv), method = adj, list = TRUE)
  posthoc.res <- data.frame(t(data.frame(posthoc$P.adjusted)))
  names(posthoc.res) <- posthoc$comparisons
    
    all.res = data.frame(test.id = asv)
    all.res$kw.p.value <- kw.result$p.value
    all.res <- cbind(all.res, posthoc.res)
    return(all.res)
    }
  
  
###############################################
###############################################
###############################################
## RUN THE DUNN 

run.dunn <- function(mat, metadata, independent.variable, adj){
  
  dunn.results <- list()
  for (i in 1:nrow(mat)){
    result <- dunn(mat[i,], metadata, independent.variable, adj)
    dunn.results[[i]] <- as.data.frame(result)
  }
  
  dunn.results <- bind_rows(dunn.results)
  
  ## now correct for multiple testing on multiple bacteria
  new.p = list()
  for(i in 2:ncol(dunn.results)){
    new.p[[i]] <- as.data.frame(p.adjust(dunn.results[,i], method = "BH"))
    
  }
  
  new.res <- bind_cols(new.p)
  final.res <- cbind(dunn.results[,1], new.res)
  names(final.res) <- names(dunn.results)
  
  return(final.res)
}

################################################
################################################

########## Friedman Test - Like a paired Kruskal wallis 
########################################################


fried <- function(row, metadata, independent.variable, sample.variable, comparisons=4){
  
  bug.name <- rownames(row)
  bug.df <- data.frame(Sample=names(row), Rel_Abun=unlist(row))
  
  metadata <- metadata[match(bug.df$Sample, rownames(metadata)), ]

  bug.df$independent.variable <- metadata[,independent.variable]
  bug.df$sample <- metadata[,sample.variable]
  
  bug.df <- bug.df[with(bug.df, order(sample, independent.variable)),]
  
  
  # Set as factors
  bug.df$independent.variable <- as.factor(bug.df$independent.variable)
  
  
  ## Friedman
  fried.res <- friedman.test(Rel_Abun ~ independent.variable | sample, data = bug.df)
  res <- data.frame(Bacteria = bug.name, fried.p = fried.res$p.value)
  
  ## Pairwise
  pairwise <- pairwise.wilcox.test(bug.df$Rel_Abun, bug.df$independent.variable, paired = TRUE, p.adjust.method = "BH")
  pairwise.res <- data.frame(melt(pairwise[[3]]))
  
  names(pairwise.res) <- c("variable.1", "variable.2", "p.value")
  # Remove the repeat or same:same tests
  pairwise.res <- pairwise.res[!is.na(pairwise.res$p.value),]
  
  pairwise.res$Bacteria <- rep(bug.name, nrow(pairwise.res))
  
  all.res <- list(res, pairwise.res)
  names(all.res) <- c("friedman", "pairwise")
  return(all.res)
  
}

##############################################################
##############################################################
#### Execute on multiple rows in a counts table ##############
##############################################################
##############################################################


run.fried <- function(normalised.counts, metadata, independent.variable="Condition", sample.variable="SampleID", comparisons=4){
  
  fried.res <- list()
  pairwise.res <- list()
  
  for (i in 1:nrow(normalised.counts)){
    result <- fried(normalised.counts[i,], metadata, independent.variable, sample.variable, comparisons=4)
    fried.res[[i]] <- as.data.frame(result$friedman)
    pairwise.res[[i]] <- as.data.frame(result$pairwise)
  }
  
  all.fried.res <- bind_rows(fried.res)
  rownames(all.fried.res) <- rownames(normalised.counts)
  all.fried.res$padj <- p.adjust(all.fried.res$fried.p)
  all.fried.res$padjstar <- stars.pval(all.fried.res$padj)
  
  # sig.bugs <- all.fried.res[all.fried.res$padj < 0.05, 1]
  
  all.pairwise.res <- bind_rows(pairwise.res)
  # all.pairwise.res <- all.pairwise.res[all.pairwise.res$Bacteria %in% sig.bugs,]
  all.pairwise.res$pstar <- stars.pval(all.pairwise.res$p.value)
  all.pairwise.res$padj <- p.adjust(all.pairwise.res$p.value, method = "BH")
  all.pairwise.res$padjstar <- stars.pval(all.pairwise.res$padj)
  
  all.res <- list(all.fried.res, all.pairwise.res)
  names(all.res) <- c("friedman", "pairwise")
 
  return(all.res) }


#####################################################
############## Wilcox test ##########################

########## Wilcox.test 
########## Comparison between 2 variables

wil.test <- function(row, metadata, independent.variable){
  iv <- unlist(metadata[, independent.variable])
  wil.result <- wilcox.test(unlist(row) ~ as.factor(iv))
  result <- data.frame(wilcox.p = wil.result$p.value)
  result$Bacteria <- rownames(row)
  return(result)}


########### Run over a data frame by row

run.wil.test <- function(df, metadata, independent.variable){
  wil.results = list()
  for(i in 1:nrow(df)){
    result <- wil.test(df[i,], metadata, independent.variable)
    wil.results[[i]] <- as.data.frame(result)
    
  }
  wil.results <- bind_rows(wil.results)
  wil.results$padj <- p.adjust(wil.results$wilcox.p, method = "BH")
  wil.results$padj.star <- stars.pval(wil.results$padj)
  return(wil.results)}



