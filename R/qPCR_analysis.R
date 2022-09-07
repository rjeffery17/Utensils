library(dplyr)



##### For inferring concentration from Ct values 
##### qPCR data with standard curve


### Concentration column must be named "Concentration"
### CT columns in std and sample data must be called "CT"

pred.conc <- function(std.data, sample.data){
  
  # Remove blank
  std.data <- std.data[!std.data$Concentration == 0,]
  std.data$log.conc <- log10(std.data$Concentration)
  lm.eq <- lm(log.conc ~ CT, data = std.data)
  
  sample.ct <- data.frame(CT = sample.data$CT)
  pred.test <- 10^(predict(lm.eq, newdata=sample.ct))
  sample.data$Concentration <- pred.test
  
  return(sample.data)
}

##### Plot Std Curve + equation for conc fitting

plot.std.curve <- function(std.data){
  # Remove blank
  std.data <- std.data[!std.data$Concentration == 0,]
  std.data$log.conc <- log10(std.data$Concentration)
  lm.eq <- lm(log.conc ~ CT, data = std.data)
  intercept <- format(lm.eq$coefficients[1], digits=4)
  Ctcoef <- format(lm.eq$coefficients[2], digits = 4)
  
  
  eq <- paste0("LC = ", "(", Ctcoef, "*Ct)", " + ", intercept)
  
  plot <- ggplot(std.data, aes(x=Concentration, y = CT)) +
    geom_point(alpha=0.5)+
    scale_x_log10() +
    theme_classic() +
    theme(text=element_text(size=15)) +
    geom_smooth(method=lm, se = FALSE) +
    geom_text(x = 0, y = 31, label = eq) +
    labs(title = std.data$Bug)
  
  return(plot)
  
}

#################################

bug.iterator <- function(bug.list, df){
  dfs <- list()
  for(i in 1:length(bug.list)){
    new_df <- df[df$Bug == bug.list[i],]
    dfs[[i]] <- new_df}
  return(dfs)
  
}

################################
gen.plot.list <- function(df.list){
  plot.list <- list()
  for(i in 1:length(df.list)){
    plot <- plot.std.curve(df.list[[i]])
    plot.list[[i]] <- plot
  }
  return(plot.list)
}
################################
pred.conc.XTREME <- function(std_df_list, sample_df_list){
  results <- list()
  for(i in 1:length(std_df_list)){
    for(i in 1:length(sample_df_list)){
      result <- pred.conc(std_df_list[[i]], sample_df_list[[i]])
      results[[i]] <- result
    }}
  results <- bind_rows(results)
  return(results)
}
  
  