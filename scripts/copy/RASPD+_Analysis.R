#!/usr/bin/env Rscript
# This script helps to filter and select compounds output by a  RASPD+ virtual screening by applying filters.
# Author(s): Jonathan Teuffel, Heidelberg Institute for Theoretical Studies (HITS) and Heidelberg University, Germany
# Licensed under the EUPL 1.2
# Requirements: 
#    Latest version of R 
#    corrplot package (https://cran.r-project.org/web/packages/corrplot/vignettes/corrplot-intro.html) 
#    RASPD+ output file: FinalResult.txt
#    Target.smi file: SMILES strings of compounds to be filtered
#start
#########Automated Script to analyze and filter RASPD+-results####################
######################Version: 12.03.2020 II######################################
##########################Jonathan Teuffel#####################

#clear workspace to remove potenital artifacts
rm(list = ls())

library(methods)

#Load or install&load all required Packages
packages = c("corrplot", "readr")
package.check <- lapply(
  packages,
  FUN = function(x) {
    if (!require(x, character.only = TRUE)) {
      install.packages(x, dependencies = TRUE)
      library(x, character.only = TRUE)
    }
  }
)
#read input-argument
args <- commandArgs(trailingOnly = T)

#move to directory
setwd("./")

#remove old outputs from the working-directory
#Define the file-names that will be deleted
old_outp <- c("scores_pairwise.png","scores_distribution.png","scores_boxplot.png","corplot_colour.png",
              "mols_pass.smi", "mols_pass.csv","Score-distribution_with_threshold.png")

rm_files <- c()
for (i in 1:length(old_outp)){
  if(file.exists(old_outp[i]) == T){
    rm_files <- old_outp[i]
    file.remove(rm_files)
  }
}
rm(old_outp,rm_files)

#check if input-files exist, if not: quit the script and give an error
if (!(file.exists("FinalResult.txt"))) {
  cat("File: 'FinalResult.txt' not found, quitting... \n")
  quit(save="no")
}
if (!(file.exists("target.smi"))) {
  cat("File: 'target.smi' not found, quitting... \n")
  quit(save="no")
}

#read FinalResult.txt
data <-  read.csv("FinalResult.txt", sep = ";", header = T)
#check if the dataset only features one score-distribution
if(dim (data)[2] == 2){
  args <- as.numeric(args)
  #If no sigma-cutoff is provided: set it to 1.5*Sigma (performed best in testing)
  if (length(args) == 0 | isTRUE(is.na(args[1]))){
    args <- c(1.5)
    print("using default-cutoff: 1.5 sigma")
  }
  #Set cutofff to x*sigma and calculate this value
  cutoff <- sqrt(var(data[,2]))
  cat("Input is only one score distribution\n cutoff is: ",args," * Sigma\n ")
  cutoff <- (args[1]*cutoff)
  cutoff <- (mean(data[,2])-cutoff)
  #Export a .png-file of the distribution with a vertical line at the custom threshold
  png(file="Score-distribution_with_threshold.png",
      width=750, height=763)
  par(cex.axis = 1.5, cex.main = 2.5, cex.lab = 1.6)
  hist(data[,2], xlab = "Predicted Binding Free Energy [kcal/mol]", ylab = "Frequency",
       main = "Score-Distribution with Cut-Off", breaks = 100)
  abline(v=cutoff,col="red", lwd = 2.5)
  legend(-7, y=20000L, legend = c("Cut-off: scores on the left\nof the line pass the filter"),
         col = "red",cex = 1.1,lty=1:0)
  dev.off()
  #Generate a new dataframe data_pass which only contains the molecules which passed the filter
  data_pass <- subset(data, data[,2] < cutoff)
  cat(dim(data_pass)[1],"molecules passed the threshold\n")
  #Read SMILES from target.smi and convert them into a data.frame
  mols = read.csv("target.smi", header = F)
  mols <- data.frame(c(seq(1:length(data[,1]))),mols)
  #create new data frame mols_pass which only contains the strings which passed the fiter
  check <- rep(F, (dim(data)[1]))
  check[data_pass$MoleculeID] <- T
  mols_pass <- subset(mols,check[mols[,1]] == T)
  #Export smiles-file and .csv file that contains the smiles-descriptor as well as MolID &score
  export_mol <- data.frame(mols_pass$V1)
  write_csv(export_mol,"mols_pass.smi", quote = F, col_names = F)
  data <- data.frame(data_pass$MoleculeID, export_mol ,data_pass[,2])
  colnames(data) <- c("ID","SMILES-string","score")
  write.csv2(data,"mols_pass.csv", quote = F)
}
#check if the dataset features seven score-distributions
if(dim(data)[2] == 8){
  # test if the arguments make sense if not, apply default: both filters
  if (isTRUE(args[1] != "f1") & isTRUE(args[1] != "f2")| is.na(args[1])) {
    print("Input is not f1 or f2 --> using defaults: applying both filters")
    args[1] <- "f2"
  }
  #If the dataset is very large: reduce scope of the dataset TO 1/100
  #for pairwise score-plots (R crashes if millions of points are plotted)
  if(isTRUE(length(data$MoleculeID) > 100000)){
  rand <- floor(runif(length(data$MoleculeID), min=0, max=101))
  data_plot <- data.frame(data,rand)
  data_plot <- subset(data, rand == 1)
  data_plot <- data_plot[,2:8]
  colnames(data_plot) <- c("ERF", "RF", "DNN", "KNN", "ISVR", "SVR", "LR")
  #export paiwise correlations-plot matrix
  png(file="scores_pairwise.png",
      width=750, height=763)
  par(cex.lab = 1.75, cex.axis = 1.75)
  plot(data_plot, cex = 0.1, main = "Pairwise Score-Correlations between All Regressions")
  dev.off()
  #clean up
  rm(data_plot)
  rm(rand)
  }
  #plot all points if the dataset is small enough to do so
  else{
  data_plot <- data[,2:8]
  colnames(data_plot) <- c("ERF", "RF", "DNN", "KNN", "ISVR", "SVR", "LR")
  #export paiwise correaltions-plot matrix
  png(file="scores_pairwise.png",
      width=750, height=763)
  par(cex.lab = 1.75, cex.axis = 1.75)
  plot(data_plot, cex = 0.1, main = "Pairwise Score-Correlations between All Regressions")
  dev.off()
  #clean up
  rm(data_plot)
  }
  #export score-distributions for each model
  png(file="scores_distribution.png",
                width=1000, height=1000)
  par(mfrow = c(2,4),cex.lab = 1.5, cex.axis = 1.5, cex.main = 3)
  hist(data$PBFE.ERF., breaks = 50, xlab = "Predicted Binding Free Energy [kcal/mol]", ylab = "Frequency", main = "ERF")
  hist(data$PBFE.RF., breaks = 50, xlab = "Predicted Binding Free Energy [kcal/mol]", ylab = "Frequency", main = "RF")
  hist(data$PBFE.DNN., breaks = 50, xlab = "Predicted Binding Free Energy [kcal/mol]", ylab = "Frequency", main = "DNN")
  hist(data$PBFE.KNN., breaks = 50, xlab = "Predicted Binding Free Energy [kcal/mol]", ylab = "Frequency", main = "KNN")
  hist(data$PBFE.SVR., breaks = 50, xlab = "Predicted Binding Free Energy [kcal/mol]", ylab = "Frequency", main = "ISVR")
  hist(data$PBFE.εSVR., breaks = 50, xlab = "Predicted Binding Free Energy [kcal/mol]", ylab = "Frequency", main = "SVR")
  hist(data$PBFE.LR., breaks = 50, xlab = "Predicted Binding Free Energy [kcal/mol]", ylab = "Frequency", main = "LRR")
  plot(c(0, 1), c(0, 1), ann = F, bty = 'n', type = 'n', xaxt = 'n', yaxt = 'n')
  text(x = 0.5, y = 0.5, paste("Frequency\n of\n predicted\n binging-free\n energies\n for each\n regression-\n model\n",
                               "i.e. score-\n distributions"), 
       cex = 2.3, col = "black")
  dev.off()
  #save pearson-coefficient-matrix for later
  res1 <- cor(data[,2:8])
  #compute quantiles for each score-distribution and set threshold to lowest quartile 
  thresholds <- c()
  for(i in 2:8){
    quart <- quantile(data[,i])
    thresholds[i] <- quart[2]
  }
  #save non-filtered data for later
  datacor1 = data[,2:8]
  #filter 1: keep only entries with only scores which are in the lowest quartile of their respective
  #score-distribution:
  data <- subset(data, data$PBFE.ERF. < thresholds[2] & data$PBFE.RF.  <  thresholds[3] & 
                   data$PBFE.DNN. < thresholds[4] & data$PBFE.KNN. <  thresholds[5] &
                   data$PBFE.SVR. < thresholds[6] & data$PBFE.εSVR. < thresholds[7] &
                   data$PBFE.LR. < thresholds[8])
  #save correlation-coefficient-matrix and scores for later
  datacor3 <- data[,2:8]
  res2 <- cor(data[,2:8])
  #filter2: Difference between scores: compute sum of the pairwise, euclid distances between scores
  if (args[1] == "f2") {
  dis2 = c()
  distance <- function(i) 
  { 
    d <- c()
    d <- abs(data[i,2]-data[i,3])+
      abs(data[i,2]-data[i,4])+
      abs(data[i,2]-data[i,5])+
      abs(data[i,2]-data[i,6])+
      abs(data[i,2]-data[i,7])+
      abs(data[i,2]-data[i,8])+
      abs(data[i,3]-data[i,4])+
      abs(data[i,3]-data[i,5])+
      abs(data[i,3]-data[i,6])+
      abs(data[i,3]-data[i,7])+
      abs(data[i,3]-data[i,8])+
      abs(data[i,4]-data[i,5])+
      abs(data[i,4]-data[i,6])+
      abs(data[i,4]-data[i,7])+
      abs(data[i,4]-data[i,8])+
      abs(data[i,5]-data[i,6])+
      abs(data[i,5]-data[i,7])+
      abs(data[i,5]-data[i,8])+
      abs(data[i,6]-data[i,7])+
      abs(data[i,6]-data[i,8])+
      abs(data[i,7]-data[i,8])
  }
  dis <- distance(c(1:length(data$MoleculeID)))
  #set threshold to lowest quartile of the distance-distribution & only keep entries with a lower distance
  cutoffs <- quantile(dis)
  cutoff <- (cutoffs[2])
  data <- data.frame(data, dis)
  data <- subset(data, dis <= (cutoff))
  cat("Using both filters\n,",dim(data)[1],"molecules passed\n")
  #save data for later
  datacor2 = data[,2:8]
  #compute extreme-values for plotting
  borders = c(min(datacor1),max(datacor1))
  #take datacor-data frames for plotting and change their col names
  colnames(datacor1) <- c("ERF", "RF", "DNN", "KNN", "ISVR", "SVR", "LR")
  colnames(datacor2) <- c("ERF", "RF", "DNN", "KNN", "ISVR", "SVR", "LR")
  #plot score-distributions before & after applying the filters
  png(file="scores_boxplot.png",
      width=750, height=750)
  par(mfrow = c(1,2),cex.lab = 1.2, cex.axis = 1, oma = c(0, 0, 2, 0))
  boxplot(datacor1, ylim = c(borders), xlab = "Scores before applying the filter(s)",ylab = "Predicted Binding Free Energy [kcal/mol]")
  boxplot(datacor2, ylim = c(borders), xlab = "Scores after applying both filters", ylab = "Predicted Binding Free Energy [kcal/mol]")
  mtext("Comparison of the Score-Distributions before and after applying the filters", outer = TRUE, cex = 1.5)
  dev.off()
  res3 <- cor(data[,2:8])
  row.names(res1) <- c("ERF", "RF", "DNN", "KNN", "ISVR", "SVR", "LR")
  colnames(res1) <- c("ERF", "RF", "DNN", "KNN", "ISVR", "SVR", "LR")
  row.names(res2) <- c("ERF", "RF", "DNN", "KNN", "ISVR", "SVR", "LR")
  colnames(res2) <- c("ERF", "RF", "DNN", "KNN", "ISVR", "SVR", "LR")
  row.names(res3) <- c("ERF", "RF", "DNN", "KNN", "ISVR", "SVR", "LR")
  colnames(res3) <- c("ERF", "RF", "DNN", "KNN", "ISVR", "SVR", "LR")
  #plot the correlation-matrices
  if(isTRUE(length(data$MoleculeID) > 1)){
  png(file="corplot_colour.png",
      width=1500, height= 650)
  par(mfrow = c(1,3),cex.lab = 1.7, cex.axis = 1.7, oma = c(0, 0, 3, 0))
  corrplot(res1, tl.col = "darkblue", tl.cex = 2, cl.cex = 2)
  corrplot(res2, tl.col = "darkblue", tl.cex = 2, cl.cex = 2)
  corrplot(res3, tl.col = "darkblue", tl.cex = 2, cl.cex = 2)
  mtext("Pearson-Correlation Matrices of the Score-predictions of the raw data (left), after applying filter f1 (middle) and after applying filter f2 (right)",
        outer = TRUE, cex = 2)
  dev.off()
  }
  #avoid plotting correlation-plots, if not enough compounds passed the filter
  else{
    cat("Skipping to generate Correlation-Plot (corplot_colour.png) due to insufficient amount of data (No. of passed compounds <= 1)\n")
  }
  }
  #Branch for when the user only wants to use the 1st filter
  if(args[1] == "f1"){
  cat("Using only the first filter,\n",dim(data)[1],"molecules passed\n")
  borders = c(min(datacor1),max(datacor1))
  colnames(datacor1) <- c("ERF", "RF", "DNN", "KNN", "ISVR", "SVR", "LR")
  colnames(datacor3) <- c("ERF", "RF", "DNN", "KNN", "ISVR", "SVR", "LR")
  png(file="scores_boxplot.png",
      width=750, height=750)
  par(mfrow = c(1,2),cex.lab = 1.2, cex.axis = 1, oma = c(0, 0, 2, 0))
  boxplot(datacor1, ylim = c(borders), xlab = "Scores before applying the filter(s)",ylab = "Predicted Binding Free Energy [kcal/mol]")
  boxplot(datacor3, ylim = c(borders), xlab = "Scores after applying the filter", ylab = "Predicted Binding Free Energy [kcal/mol]")
  mtext("Comparison of the Score-Distributions before and after applying the filter", outer = TRUE, cex = 1.5)
  dev.off()
  row.names(res1) <- c("ERF", "RF", "DNN", "KNN", "ISVR", "SVR", "LR")
  colnames(res1) <- c("ERF", "RF", "DNN", "KNN", "ISVR", "SVR", "LR")
  row.names(res2) <- c("ERF", "RF", "DNN", "KNN", "ISVR", "SVR", "LR")
  colnames(res2) <- c("ERF", "RF", "DNN", "KNN", "ISVR", "SVR", "LR")
  if (isTRUE(length(data$MoleculeID) > 1)) {
  png(file="corplot_colour.png",
      width=1500, height= 650)
  par(mfrow = c(1,2),cex.lab = 1.7, cex.axis = 1.7, oma = c(0, 0, 2, 0))
  corrplot(res1, tl.col = "darkblue", tl.cex = 2, cl.cex = 2)
  corrplot(res2, tl.col = "darkblue", tl.cex = 2, cl.cex = 2)
  mtext("Pearson-Correlation Matrices of the Score-predictions of the raw data (left) and after applying filter f1 (right)",
        outer = TRUE, cex = 2)
  dev.off()
  }
  else{
    cat("Skipping Correlation-Plots due to insufficient amount of data (No. of passed compounds <= 1)\n")
  }
  }
  #Prepare Export (see: first part of the script)
  if (isTRUE(length(data$MoleculeID) > 1)){
    mols = read.csv("target.smi", header = F)
    mols <- data.frame(c(seq(1:(dim(datacor1)[1]))),mols)
    colnames(mols) <- c("ID", "SMILE")
    check <- rep(F, (dim(data)[1]))
    check[data$MoleculeID] <- T
    mols_pass <- subset(mols, check[mols$ID] == T)
    #Export-part (see: first part of the script)
    export_mol <- data.frame(mols_pass$SMILE)
    write_csv(export_mol,"mols_pass.smi", quote = F, col_names = F)
    data <- data.frame(data$MoleculeID, export_mol ,data$PBFE.ERF., data$PBFE.RF., data$PBFE.DNN. ,
                       data$PBFE.KNN., data$PBFE.SVR., data$PBFE.εSVR., data$PBFE.LR.)
    colnames(data) <- c("ID","SMILE-descriptor","ERF","RF","DNN","KNN","ISVR","εSV","LR")
    write.csv2(data,"mols_pass.csv", quote = F)
  }
  #avoid errors when trying to create a data-frame from only 1 result
  else{
    cat("Insufficient number of passed compounds for data-export \n")
  }
}
