#Nov 21 2019. lmer from lmerTest packge for analysis of relative abundances (log10(1+1E-07))

library(lmerTest)
setwd("~/Dropbox/DPhil/CD1 Maternal Transfer/16S sequencing/CMS008")

#read in untransformed relative abundance data; there are 7 columns with sample info before the data starts
rel_abun <- read.csv("samples_otus_rel.csv")
#log transform
rel_abun_log <- data.frame(lapply(rel_abun[8:length(rel_abun)],function(x) log10(x + 1e-07)))
rel_abun_log <- cbind(rel_abun[1:7], rel_abun_log)
#subset to offspring only
rel_abun_log <- subset(rel_abun_log, Cohort == "F1")
#remove extra sample info columns; remaining data frame has first column as "Group" and second column as "Litter"
rel_abun_log <- rel_abun_log[,-c(1:3,5:6)]

#my very inelegant way of doing lmer on each ASV
#significance matrix
sig <- data.frame(matrix(nrow = ncol(rel_abun_log)-2, ncol = 1))
rownames(sig) <- colnames(rel_abun_log[3:length(rel_abun_log)]) 
colnames(sig) <- c("P")

#matrix for lmer for each ASV
#make dataframe that the lmer test will be run on
df <- data.frame(matrix(ncol = 1, nrow = nrow(rel_abun_log)))
df <- cbind(rel_abun_log$Group, rel_abun_log$Litter, df)
names(df) <- c("Group","Litter","Abundance")

for (i in 3:length(rel_abun_log)) {
  df[3] <- rel_abun_log[i]
  names(df)[3] <- c("Abundance")
  lmer <- lmer(Abundance ~ Group + (1|Litter), data = df)
  summary <- anova(lmer)
  sig[i-2,1] <- summary$P
  }

library(dplyr)
library(tibble)
sig_sorted <- sig
sig_sorted <- rownames_to_column(sig_sorted)
sig_sorted <- arrange(sig_sorted,P)
