#!/usr/bin/env Rscript

#Load required packages

rm(list = ls())
args = commandArgs(trailingOnly=TRUE)

library(ggplot2)
library(dplyr)
library(DescTools)
library(pROC)
library(caret)

##############################
#######Read in the files

out_PRSice = args[1]
out_lasso = args[2]
out_PRScs = args[3]
out_LDpred2 = args[4]
out_lasso2 = args[5]
out_comparison = args[6]

test_prefix = args[7]
test_prefix_rsID = args[8]
training_prefix = args[9]
training_prefix_rsID = args[10]

cov_test = args[11]
cov_training = args[12]
pheno_file = args[13]

lasso_thresholds = args[14]
phi_values = args[15]

cross_validation = args[16]

##############################
#######Define functions

#### PC correction

PC_correction <- function(data, score, score_res) {
  tip <- paste(score, "~", PCs, sep="")
  alpha <- lm(tip, data = data)
  beta <- residuals(alpha, "response")
  data[[score_res]] <- beta
  return(data)
}

#### Standardize

standardize <- function(data_ref,data_target, score, score_sc) {
  mean_ref <- mean(data_ref[[score]])
  sd_ref <- sd(data_ref[[score]])
  data_target[[score_sc]] <- (data_target[[score]] - mean_ref)/sd_ref
  return(data_target)
}

#### Figures

##### Barplot_gradient

barplot_gradient <- function(data, Pt, R2, P_value, all_scores, fill_color, out, Tool){
  if(length(unique(data[[P_value]])) == 1){
    bar_plot <- ggplot(data=data, aes(x=!!sym(Pt), y=!!sym(R2), fill=fill_color)) +
      geom_bar(stat = "identity", position = "dodge") +
      scale_x_discrete(limits=all_scores, guide = guide_axis(angle = 90)) +
      scale_fill_identity(name="-log10(P_value)", guide="none") +
      xlab("") + ylab("R2") + theme_minimal() +
      ggtitle(Tool)
  } else {
    bar_plot <- ggplot(data=data, aes(x=!!sym(Pt), y=!!sym(R2), fill=-log10(P_value))) +
      geom_bar(stat = "identity", position = "dodge") +
      scale_x_discrete(limits=all_scores, guide = guide_axis(angle = 90)) +
      scale_fill_gradient2(high = fill_color, name="-log10(P_value)") +
      xlab("") + ylab("R2") + theme_minimal() +
      ggtitle(Tool)
  }
  ggsave(out, bar_plot, width=7, height=7)
}


#### Correlation plot

correlation_plot <- function(data, x_var, y_var, point_color, out){
  #Compute correlation
  spearman_cor <- cor(data[[x_var]], data[[y_var]], method = "spearman")
  spearman_r2 <- round(spearman_cor^2, 3)
  #Determine position for the annotation
  x_min <- min(data[[x_var]], na.rm = TRUE)
  y_max <- max(data[[y_var]], na.rm = TRUE)
  
  #Create the plot
  plot <- ggplot(data=data, aes(x=!!sym(x_var), y=!!sym(y_var))) +
    geom_point(color=point_color) +
    geom_smooth(method="lm", se=F, color="grey", formula = y~x) +
    annotate("text", x=x_min, y=y_max, label=paste0("Spearman R² = ", spearman_r2),
             hjust=0, vjust=1, size=4, color="black") +
    theme_minimal() +
    xlab("PRS") + ylab("PHENO")
  
  #Save the plot
  ggsave(filename=out, plot=plot, height=6, width=9)
}


##################
##General files###
##################

PCs_test <- read.table(cov_test, header=T)
PCs_training <- read.table(cov_training, header=T)

PC_only <- select(PCs_test, -FID, -IID)
PCs <- paste(colnames(PC_only), collapse = '+')

pheno <- read.table(pheno_file, header=T)
colnames(pheno) <- c("FID", "IID", "PHENO")
pheno <- pheno[!is.na(pheno$PHENO),]
pheno$PHENO <- as.numeric(pheno$PHENO)


##################
#####PRSice#######
##################

cat("Start get best PRS... PRSice-2...\n")

#### Scores

PRS_test_PRSice <- read.table(paste0(out_PRSice, "/", test_prefix, ".all_score"), header=T)
PRS_training_PRSice <- read.table(paste0(out_PRSice, "/", training_prefix, ".all_score"), header=T)

##############################
#######PC correction

#### Get Pt

get_Pt <- select(PRS_test_PRSice, -FID, -IID)
Pt <- colnames(get_Pt)
Pt_res <- paste(Pt, "_res", sep="")

#### Merge files with PCs

PRS_test_PRSice_PC <- merge(PRS_test_PRSice, PCs_test, by=c("IID", "FID"))
PRS_training_PRSice_PC <- merge(PRS_training_PRSice, PCs_training, by=c("IID", "FID"))

#### Merge test and training for PC correction

PRS_test_PRSice_PC$Group <- "test"
PRS_training_PRSice_PC$Group <- "training"
merged_PRSice <- rbind(PRS_test_PRSice_PC, PRS_training_PRSice_PC)
##### Account for NA columns
na_col_PRSice <- colSums(is.na(merged_PRSice))
all_na_cols_PRSice <- which(na_col_PRSice == nrow(merged_PRSice))
if(length(all_na_cols_PRSice) > 0) {
  merged_PRSice <- merged_PRSice[, -all_na_cols_PRSice]
} else {
  # If all_na_cols is empty, just use the original dataframe
  merged_PRSice <- merged_PRSice
}
##### Account for rows that have all NA
na_rows_PRSice <- apply(merged_PRSice, 1, function(row) any(is.na(row)))
merged_PRSice <- merged_PRSice[!na_rows_PRSice,]

#### PC correction

for (i in 1:length(Pt)) {
  score <- Pt[i]
  score_res <- Pt_res[i]
  merged_PRSice <- PC_correction(data=merged_PRSice, score=score, score_res=score_res)
}

#### Split again to then do standardization
PRS_test_PRSice_PC <- merged_PRSice[(merged_PRSice$Group=="test"),]
PRS_test_PRSice_PC <- PRS_test_PRSice_PC[, !colnames(PRS_test_PRSice_PC) %in% "Group"]
PRS_training_PRSice_PC <- merged_PRSice[(merged_PRSice$Group=="training"),]
PRS_training_PRSice_PC <- PRS_training_PRSice_PC[, !colnames(PRS_training_PRSice_PC) %in% "Group"]


##############################
#######Standardization

#### Standardize

Pt_sc <- c(Pt, Pt_res)
Pt_sc_2 <- paste(Pt_sc, "_sc", sep="")

for (i in 1:length(Pt_sc)) {
  score <- Pt_sc[i]
  score_sc <- Pt_sc_2[i]
  PRS_test_PRSice_PC <- standardize(data_ref=PRS_training_PRSice_PC, data_target=PRS_test_PRSice_PC, score=score, score_sc=score_sc)
  PRS_training_PRSice_PC <- standardize(data_ref=PRS_training_PRSice_PC, data_target=PRS_training_PRSice_PC, score=score, score_sc=score_sc)
}

#### Get only scores to write to a file

PRS_test_PRSice_scores <- colnames(PRS_test_PRSice_PC)[!(colnames(PRS_test_PRSice_PC) %in% colnames(PC_only))]
PRS_test_PRSice <- PRS_test_PRSice_PC[, PRS_test_PRSice_scores]
PRS_training_PRSice_scores <- colnames(PRS_training_PRSice_PC)[!(colnames(PRS_training_PRSice_PC) %in% colnames(PC_only))]
PRS_training_PRSice <- PRS_training_PRSice_PC[, PRS_training_PRSice_scores]

write.table(PRS_test_PRSice, paste0(out_PRSice, "/", test_prefix, "_scaled_scores"), row.names = F, quote = F, sep="\t")
write.table(PRS_training_PRSice, paste0(out_PRSice, "/", training_prefix, "_scaled_scores"), row.names = F, quote = F, sep="\t")


##############################
#######Linear regression

Pt_res_sc <- paste(Pt, "_res_sc", sep="")
PRS_test_PRSice_pheno <- merge(PRS_test_PRSice, pheno, by=c("IID", "FID"))
PRS_test_PRSice_pheno$PHENO <- as.numeric(PRS_test_PRSice_pheno$PHENO)

#### Perform the regression

if (cross_validation=="FALSE"){
  Regression_results_PRSice <- data.frame(matrix(ncol=6, nrow=0))
  
  for (i in 1:length(Pt_res_sc)) {
    a <- Pt_res_sc[i]
    thresh <- Pt[i]  
    tip <- paste("PHENO", "~", a,sep="")
    alpha <- lm(tip, data = PRS_test_PRSice_pheno)
    R2 <- summary(alpha)$adj.r.squared
    beta <- coef(summary(alpha))[2,1]
    SE <- coef(summary(alpha))[2,2]
    p <- coef(summary(alpha))[2,4]
    Regression_results_PRSice <- rbind(Regression_results_PRSice, c("PRSice",thresh,beta,SE,p,R2))
    colnames(Regression_results_PRSice) <- c("Tool","Parameters", "Beta", "SE", "P_value", "R2")
  }
  
  Regression_results_PRSice$R2 <- as.numeric(Regression_results_PRSice$R2)
  Regression_results_PRSice$P_value <- as.numeric(Regression_results_PRSice$P_value)
  Regression_results_PRSice$P_value <- ifelse(Regression_results_PRSice$P_value==0, 3.4e-314, Regression_results_PRSice$P_value)
  
  write.table(Regression_results_PRSice, paste0(out_PRSice, "/Regression_results_PRSice"), row.names=F, quote=F, sep="\t")
  
} else if (cross_validation=="TRUE"){
  Regression_results_PRSice <- data.frame(matrix(ncol=7, nrow=0))
  train_control <- trainControl(method = "cv", number = 10)
  
  for (i in 1:length(Pt_res_sc)) {
    a <- Pt_res_sc[i]
    thresh <- Pt[i]  
    tip <- paste("PHENO", "~", a,sep="")
    alpha <- lm(tip, data = PRS_test_PRSice_pheno)
    R2 <- summary(alpha)$adj.r.squared
    beta <- coef(summary(alpha))[2,1]
    SE <- coef(summary(alpha))[2,2]
    p <- coef(summary(alpha))[2,4]
    model <- train(as.formula(tip), data=PRS_test_PRSice_pheno, method="lm", trControl=train_control, metric="Rsquared")
    R2_CV <- (model$results$Rsquared)
    Regression_results_PRSice <- rbind(Regression_results_PRSice, c("PRSice",thresh,beta,SE,p,R2,R2_CV))
    colnames(Regression_results_PRSice) <- c("Tool","Parameters", "Beta", "SE", "P_value", "R2","R2_CV")
  }
  
  Regression_results_PRSice$R2 <- as.numeric(Regression_results_PRSice$R2)
  Regression_results_PRSice$R2_CV <- as.numeric(Regression_results_PRSice$R2_CV)
  Regression_results_PRSice$P_value <- as.numeric(Regression_results_PRSice$P_value)
  Regression_results_PRSice$P_value <- ifelse(Regression_results_PRSice$P_value==0, 3.4e-314, Regression_results_PRSice$P_value)
  
  write.table(Regression_results_PRSice, paste0(out_PRSice, "/Regression_results_PRSice"), row.names=F, quote=F, sep="\t")
  
} else {
  "Fill in TRUE or FALSE for cross-validation"
}

#### Make bar plot

if (cross_validation=="FALSE"){
  barplot_gradient(data=Regression_results_PRSice, Pt="Parameters", R2="R2", P_value="P_value", all_scores=Pt, fill_color="cadetblue3", out=paste0(out_PRSice, "/bar_plot_R2_PRSice.svg"), Tool="PRSice")
} else if (cross_validation=="TRUE"){
  barplot_gradient(data=Regression_results_PRSice, Pt="Parameters", R2="R2_CV", P_value="P_value", all_scores=Pt, fill_color="cadetblue3", out=paste0(out_PRSice, "/bar_plot_R2_PRSice.svg"), Tool="PRSice")
} else {
  "Fill in TRUE or FALSE for cross-validation"
}

 
##############################
#######Select best PRS and plot

#### Get best PRS
Regression_results_valid_PRSice <- Regression_results_PRSice[(Regression_results_PRSice$P_value < 0.05 & Regression_results_PRSice$Beta > 0),]

if (cross_validation=="FALSE"){
  ##### Check if Regression_results_valid is empty
  if (nrow(Regression_results_valid_PRSice) == 0) {
    sorted_regression_results_PRSice <- data.frame(Tool = "PRSice",Parameters = colnames(PRS_test_PRSice_pheno[3]), Beta = 0, SE = 0, P_value = 1, R2 = 0.001)
    cat("There are no regression results with p < 0.05 and beta > 0... creating dummy row")
  } else { 
    sorted_regression_results_PRSice <- Regression_results_valid_PRSice %>% arrange(desc(R2))
  }
  
  best_Pt <- sorted_regression_results_PRSice[1,2]
  best_Pt <- paste0(best_Pt, "_res_sc")
  col_select_PRSice <- c("FID","IID", best_Pt, "PHENO")
  best_PRS_PRSice <- select(PRS_test_PRSice_pheno, all_of(col_select_PRSice))
  best_PRS_PRSice$Tool <- "PRSice"
  best_PRS_PRSice$parameters <- sorted_regression_results_PRSice[1,2]
  colnames(best_PRS_PRSice) <- c("FID","IID", "score", "PHENO", "Tool", "parameters")
  
  Regression_best_per_tool <- data.frame()
  Regression_best_PRSice <- sorted_regression_results_PRSice[1,]
  Regression_best_per_tool <- rbind(Regression_best_per_tool, Regression_best_PRSice)
  
  #### Write to new file
  
  col_select_PRSice <- c("FID","IID", "score", "Tool", "parameters")
  selection_best_PRS_PRSice <- select(best_PRS_PRSice, all_of(col_select_PRSice))
  write.table(selection_best_PRS_PRSice, paste0(out_PRSice, "/best_PRS_PRSice"), row.names=F, quote=F, sep="\t")
  
} else if (cross_validation=="TRUE"){
  ##### Check if Regression_results_valid is empty
  if (nrow(Regression_results_valid_PRSice) == 0) {
    sorted_regression_results_PRSice <- data.frame(Tool = "PRSice",Parameters = colnames(PRS_test_PRSice_pheno[3]), Beta = 0, SE = 0, P_value = 1, R2 = 0.001, R2_CV=0.001)
    cat("There are no regression results with p < 0.05 and beta > 0... creating dummy row")
  } else { 
    sorted_regression_results_PRSice <- Regression_results_valid_PRSice %>% arrange(desc(R2_CV))
  }
  
  best_Pt <- sorted_regression_results_PRSice[1,2]
  best_Pt <- paste0(best_Pt, "_res_sc")
  col_select_PRSice <- c("FID","IID", best_Pt, "PHENO")
  best_PRS_PRSice <- select(PRS_test_PRSice_pheno, all_of(col_select_PRSice))
  best_PRS_PRSice$Tool <- "PRSice"
  best_PRS_PRSice$parameters <- sorted_regression_results_PRSice[1,2]
  colnames(best_PRS_PRSice) <- c("FID","IID", "score", "PHENO", "Tool", "parameters")
  
  Regression_best_per_tool <- data.frame()
  Regression_best_PRSice <- sorted_regression_results_PRSice[1,]
  Regression_best_per_tool <- rbind(Regression_best_per_tool, Regression_best_PRSice)
  
  #### Write to new file
  
  col_select_PRSice <- c("FID","IID", "score", "Tool", "parameters")
  selection_best_PRS_PRSice <- select(best_PRS_PRSice, all_of(col_select_PRSice))
  write.table(selection_best_PRS_PRSice, paste0(out_PRSice, "/best_PRS_PRSice"), row.names=F, quote=F, sep="\t")
  
} else {
  "Fill in TRUE or FALSE for cross-validation"
}


#### Correlation best score and phenotype

correlation_plot(data=PRS_test_PRSice_pheno, x_var=best_Pt, y_var="PHENO", point_color="mediumpurple", out=paste0(out_PRSice, "/Correlation_plot_PRSice.svg"))

cat("Done getting best PRS for PRSice-2.... Lassosum....\n")

##################
#####Lassosum#####
##################

##############################
#######Read in the files

thresholds <- strsplit(lasso_thresholds, ",")
shrinkage <- character(length = length(thresholds[[1]]))
for (i in 1:length(thresholds[[1]])) {
  shrinkage[i] <- paste("s_",thresholds[[1]][i], sep = "")
}


data_frames_training <- list()
data_frames_test <- list()
for (i in 1:length(shrinkage)) {
  a <- shrinkage[i]
  file_name <- paste(out_lasso, "/", test_prefix, "_scores_", a, ".txt", sep="")
  data <- read.table(file_name, header=T)
  data <- data[, colSums(data==0) !=(nrow(data))] #Delete columns where all values are 0
  data_frames_test[[i]] <- data
  
  file_name <- paste(out_lasso, "/", training_prefix, "_scores_", a, ".txt", sep="")
  data <- read.table(file_name, header=T)
  data <- data[, colSums(data==0) !=(nrow(data))]
  data_frames_training[[i]] <- data
}

##############################
#######PC correction

data_frames_training_PC <- list()
data_frames_test_PC <- list()
PC_names <- colnames(PC_only)

for (j in 1:length(data_frames_training)){
  training <- data_frames_training[[j]]
  training_PC <- merge(training, PCs_training, by=c("FID","IID"))
  training_PC$Group <- "training"
  test <- data_frames_test[[j]]
  test_PC <- merge(test, PCs_test, by=c("FID","IID"))
  test_PC$Group <- "test"
  common_cols <- intersect(names(training_PC), names(test_PC))
  training_PC <- training_PC[, common_cols]
  test_PC <- test_PC[, common_cols]
  merge_training_test <- rbind(training_PC, test_PC)
  na_col_lasso <- colSums(is.na(merge_training_test))
  all_na_cols_lasso <- which(na_col_lasso == nrow(merge_training_test))
  if(length(all_na_cols_lasso) > 0) {
    merge_training_test <- merge_training_test[, -all_na_cols_lasso]
  } else {
    # If all_na_cols is empty, just use the original dataframe
    merge_training_test <- merge_training_test
  }
  na_rows_lasso <- apply(merge_training_test, 1, function(row) any(is.na(row)))
  merge_training_test <- merge_training_test[!na_rows_lasso, ]
  scores <- select(merge_training_test, -IID, -FID, -Group, -all_of(PC_names))
  lambda <- colnames(scores)
  lambda_res <- paste(lambda, "_res", sep="")
for (i in 1:length(lambda)) {
  score <- lambda[i]
  score_res <- lambda_res[i]
  merge_training_test <- PC_correction(data=merge_training_test, score=score, score_res=score_res)
}
  training_residuals <- merge_training_test[(merge_training_test$Group=="training"),]
  training_residuals <- training_residuals[, !colnames(training_residuals) %in% "Group"]
  test_residuals <- merge_training_test[(merge_training_test$Group=="test"),]
  test_residuals <- test_residuals[, !colnames(test_residuals) %in% "Group"]

  data_frames_training_PC[[j]] <- training_residuals
  data_frames_test_PC[[j]] <- test_residuals
}

##############################
#######Standardization

for (i in 1:length(data_frames_training_PC)) {
  training_PC <- data_frames_training_PC[[i]]
  test_PC <- data_frames_test_PC[[i]]
  training <- data_frames_training[[i]]
  scores <- select(training_PC, -IID, -FID, -all_of(PC_names))
  lambda_sc <- colnames(scores)
  lambda_sc_2 <- paste(lambda_sc, "_sc", sep="")
  
  for (j in 1:length(lambda_sc)){
    score <- lambda_sc[j]
	  score_sc <- lambda_sc_2[j]
	  training_PC <- standardize(data_ref=training_PC, data_target=training_PC, score=score, score_sc=score_sc)
	  test_PC <- standardize(data_ref=training_PC, data_target=test_PC, score=score, score_sc=score_sc)
  }
  
  training_col <- colnames(training_PC)[!(colnames(training_PC) %in% colnames(PC_only))]
  training_PC <- training_PC[, training_col]
  test_col <- colnames(test_PC)[!(colnames(test_PC) %in% colnames(PC_only))]
  test_PC <- test_PC[, test_col]
  data_frames_training[[i]] <- training_PC
  data_frames_test[[i]] <- test_PC
  
  write.table(training_PC, paste0(out_lasso, "/", training_prefix, "_", shrinkage[i], "_scaled_scores"), row.names=F, quote=F, sep="\t")
  write.table(test_PC, paste0(out_lasso, "/", test_prefix, "_", shrinkage[i], "_scaled_scores"), row.names=F, quote=F, sep="\t")
}


##############################
#######Linear regression

#### Perform the regression

if (cross_validation=="FALSE"){
  Regression_results_lasso <- data.frame(matrix(ncol=8, nrow=0))
  
  for (j in 1:length(data_frames_test)){
    test <- data_frames_test[[j]]
    scores <- colnames(test)[grep("_res_sc$", colnames(test))]
    columns <- c("IID", "FID", scores)
    test_pheno <- test[, columns]
    test_pheno <- merge(test_pheno, pheno, by=c("IID", "FID"))
    colnames(test_pheno) <- c("IID", "FID",scores, "PHENO")
    test_pheno$PHENO <- as.numeric(test_pheno$PHENO)
    
    for (i in 1:length(scores)) {
      b <- scores[i]
      thresh <- paste(b)  
      tip <- paste("PHENO", "~", b,sep="")
      alpha <- lm(tip, data = test_pheno)
      R2 <- summary(alpha)$adj.r.squared
      beta <- coef(summary(alpha))[2,1]
      SE <- coef(summary(alpha))[2,2]
      p <- coef(summary(alpha))[2,4]
      Regression_results_lasso <- rbind(Regression_results_lasso, c(shrinkage[j], thresh,beta,SE,p,R2, "lassosum", paste0(shrinkage[j], "_", thresh)))
      colnames(Regression_results_lasso) <- c("Shrink", "Lambda", "Beta", "SE", "P_value", "R2", "Tool", "Parameters")
    }
  }
  
  Regression_results_lasso$R2 <- as.numeric(Regression_results_lasso$R2)
  Regression_results_lasso$P_value <- as.numeric(Regression_results_lasso$P_value)
  Regression_results_lasso$P_value <- ifelse(Regression_results_lasso$P_value==0, 3.4e-314, Regression_results_lasso$P_value)
  
  write.table(Regression_results_lasso, paste0(out_lasso, "/Regression_results_lasso_all"), row.names=F, quote=F, sep="\t")
  
} else if (cross_validation=="TRUE"){
  Regression_results_lasso <- data.frame(matrix(ncol=9, nrow=0))
  train_control <- trainControl(method = "cv", number = 10)
  
  for (j in 1:length(data_frames_test)){
    test <- data_frames_test[[j]]
    scores <- colnames(test)[grep("_res_sc$", colnames(test))]
    columns <- c("IID", "FID", scores)
    test_pheno <- test[, columns]
    test_pheno <- merge(test_pheno, pheno, by=c("IID", "FID"))
    colnames(test_pheno) <- c("IID", "FID",scores, "PHENO")
    test_pheno$PHENO <- as.numeric(test_pheno$PHENO)
    
    for (i in 1:length(scores)) {
      b <- scores[i]
      thresh <- paste(b)  
      tip <- paste("PHENO", "~", b,sep="")
      alpha <- lm(tip, data = test_pheno)
      R2 <- summary(alpha)$adj.r.squared
      beta <- coef(summary(alpha))[2,1]
      SE <- coef(summary(alpha))[2,2]
      p <- coef(summary(alpha))[2,4]
      model <- train(as.formula(tip), data=test_pheno, method="lm", trControl=train_control, metric="Rsquared")
      R2_CV <- (model$results$Rsquared)
      Regression_results_lasso <- rbind(Regression_results_lasso, c(shrinkage[j], thresh,beta,SE,p,R2, R2_CV,"lassosum", paste0(shrinkage[j], "_", thresh)))
      colnames(Regression_results_lasso) <- c("Shrink", "Lambda", "Beta", "SE", "P_value", "R2", "R2_CV","Tool", "Parameters")
    }
  }
  
  Regression_results_lasso$R2 <- as.numeric(Regression_results_lasso$R2)
  Regression_results_lasso$R2_CV <- as.numeric(Regression_results_lasso$R2_CV)
  Regression_results_lasso$P_value <- as.numeric(Regression_results_lasso$P_value)
  Regression_results_lasso$P_value <- ifelse(Regression_results_lasso$P_value==0, 3.4e-314, Regression_results_lasso$P_value)
  
  write.table(Regression_results_lasso, paste0(out_lasso, "/Regression_results_lasso_all"), row.names=F, quote=F, sep="\t")
  
} else {
  "Fill in TRUE or FALSE for cross-validation"
}


##############################
#######Select best PRS and plot

#### Get best PRS per shrinkage

Regression_results_valid_lasso <- Regression_results_lasso[(Regression_results_lasso$P_value < 0.05 & Regression_results_lasso$Beta > 0),]

if (cross_validation=="FALSE"){
  best_models_lasso <- data.frame(matrix(ncol=8, nrow=0))
  for (i in 1:length(shrinkage)){
    a <- shrinkage[i]
    data <- Regression_results_lasso[(Regression_results_lasso$Shrink==a),]
    write.table(data, paste0(out_lasso, "/Regression_results_", a), row.names = F, quote = F, sep="\t")
    data$R2 <- as.numeric(data$R2)
    barplot_gradient(data=data, Pt="Lambda", R2="R2", P_value="P_value", all_scores=scores, fill_color="cadetblue3", out=paste0(out_lasso, "/bar_plot_R2_lassosum_", a,".svg"), Tool=paste0("Lassosum: shrinkage ", a))
    data_valid <- data[(data$P_value < 0.05 & data$Beta > 0),]
    ##### Check if Regression_results_valid is empty
    if (nrow(data_valid) == 0) {
      data_sorted <- data.frame(Shrink = a, Lambda = paste0(colnames(data_frames_test[[i]])[3], "_res_sc"), Beta = 0, SE = 0, P_value = 1, R2 = 0.001, Tool = "lassosum", Parameters = paste0(a, "_", colnames(data_frames_test[[i]])[3], "_res_sc"))
      cat("There are no regression results with p < 0.05 and beta > 0 for shrinkage ", a, ".... creating dummy row")
    } else { 
      data_sorted <- data_valid %>% arrange(desc(R2))
    }
    best_models_lasso <- rbind(best_models_lasso, data_sorted[1,])
  }
  
  #### barplot with best lambda per shrinkage
  
  barplot_gradient(data=best_models_lasso, Pt="Shrink", R2="R2", P_value="P_value", all_scores=shrinkage, fill_color="cadetblue3", out=paste0(out_lasso, "/bar_plot_R2_lassosum_best_per_s.svg"), Tool="Lassosum: best per shrinkage")
  
  #### Get best PRS
  sorted_regression_results_lasso <- best_models_lasso %>% arrange(desc(R2))
  best_s <- sorted_regression_results_lasso[1,1]
  best_lambda <- sorted_regression_results_lasso[1,2]
  
  best_PRS_lasso <- read.table(paste0(out_lasso, "/", test_prefix, "_", best_s, "_scaled_scores"), header=T)
  best_PRS_lasso <- merge(best_PRS_lasso, pheno, by=c("IID", "FID"))
  col_select_lasso <- c("FID","IID", best_lambda, "PHENO")
  best_PRS_lasso <- select(best_PRS_lasso, all_of(col_select_lasso))
  best_PRS_lasso$Tool <- "lassosum"
  best_PRS_lasso$parameters <- paste0(best_s, "_", best_lambda)
  colnames(best_PRS_lasso) <- c("FID","IID", "score", "PHENO", "Tool", "parameters")
  
  write.table(best_PRS_lasso, paste0(out_lasso, "/best_PRS_lasso"), quote = F, row.names = F, sep="\t")
  
  Regression_best_lasso <- sorted_regression_results_lasso[1,]
  Regression_best_lasso <- Regression_best_lasso[1,c(3:8)]
  Regression_best_per_tool <- rbind(Regression_best_per_tool, Regression_best_lasso)
  
} else if (cross_validation=="TRUE"){
  best_models_lasso <- data.frame(matrix(ncol=9, nrow=0))
  for (i in 1:length(shrinkage)){
    a <- shrinkage[i]
    data <- Regression_results_lasso[(Regression_results_lasso$Shrink==a),]
    write.table(data, paste0(out_lasso, "/Regression_results_", a), row.names = F, quote = F, sep="\t")
    data$R2_CV <- as.numeric(data$R2_CV)
    barplot_gradient(data=data, Pt="Lambda", R2="R2_CV", P_value="P_value", all_scores=scores, fill_color="cadetblue3", out=paste0(out_lasso, "/bar_plot_R2_lassosum_", a,".svg"), Tool=paste0("Lassosum: shrinkage ", a))
    data_valid <- data[(data$P_value < 0.05 & data$Beta > 0),]
    ##### Check if Regression_results_valid is empty
    if (nrow(data_valid) == 0) {
      data_sorted <- data.frame(Shrink = a, Lambda = paste0(colnames(data_frames_test[[i]])[3], "_res_sc"), Beta = 0, SE = 0, P_value = 1, R2 = 0.001, R2_CV=0.001,Tool = "lassosum", Parameters = paste0(a, "_", colnames(data_frames_test[[i]])[3], "_res_sc"))
      cat("There are no regression results with p < 0.05 and beta > 0 for shrinkage ", a, ".... creating dummy row")
    } else { 
      data_sorted <- data_valid %>% arrange(desc(R2_CV))
    }
    best_models_lasso <- rbind(best_models_lasso, data_sorted[1,])
  }
  
  #### barplot with best lambda per shrinkage
  
  barplot_gradient(data=best_models_lasso, Pt="Shrink", R2="R2_CV", P_value="P_value", all_scores=shrinkage, fill_color="cadetblue3", out=paste0(out_lasso, "/bar_plot_R2_lassosum_best_per_s.svg"), Tool="Lassosum: best per shrinkage")
  
  #### Get best PRS
  sorted_regression_results_lasso <- best_models_lasso %>% arrange(desc(R2_CV))
  best_s <- sorted_regression_results_lasso[1,1]
  best_lambda <- sorted_regression_results_lasso[1,2]
  
  best_PRS_lasso <- read.table(paste0(out_lasso, "/", test_prefix, "_", best_s, "_scaled_scores"), header=T)
  best_PRS_lasso <- merge(best_PRS_lasso, pheno, by=c("IID", "FID"))
  col_select_lasso <- c("FID","IID", best_lambda, "PHENO")
  best_PRS_lasso <- select(best_PRS_lasso, all_of(col_select_lasso))
  best_PRS_lasso$Tool <- "lassosum"
  best_PRS_lasso$parameters <- paste0(best_s, "_", best_lambda)
  colnames(best_PRS_lasso) <- c("FID","IID", "score", "PHENO", "Tool", "parameters")
  
  write.table(best_PRS_lasso, paste0(out_lasso, "/best_PRS_lasso"), quote = F, row.names = F, sep="\t")
  
  Regression_best_lasso <- sorted_regression_results_lasso[1,]
  Regression_best_lasso <- Regression_best_lasso[1,c(3:9)]
  Regression_best_per_tool <- rbind(Regression_best_per_tool, Regression_best_lasso)
} else {
  "Fill in TRUE or FALSE for cross-validation"
}

correlation_plot(data=best_PRS_lasso, x_var="score", y_var="PHENO", point_color="mediumpurple", out=paste0(out_lasso, "/Correlation_plot_lassosum.svg"))

cat("Done getting best PRS for lassosum.... PRS-CS....\n")

##################
#####PRS-CS#######
##################

##############################
#######Read in the files

phi <- strsplit(phi_values, ",")

phi_values_raw <- character(length = length(phi[[1]]))
for (i in 1:length(phi[[1]])) {
  phi_values_raw[i] <- paste("phi_",phi[[1]][i], sep = "")
  phi_values[i] <- gsub("[+-]", ".", phi_values_raw[i])
}

PRS_test_PRS_CS <- PCs_test[,c("IID", "FID")]
PRS_training_PRS_CS <- PCs_training[,c("IID","FID")]

for (i in 1:length(phi_values_raw)) {
  a <- phi_values_raw[i]
  b <- phi_values[i]
  
  #test
  file_name <- paste(out_PRScs, "/", test_prefix, "_", a, ".profile", sep="")
  data <- read.table(file_name, header=T)
  data <- data[,c(1,2,6)]
  colnames(data) <- c("FID", "IID", b)
  PRS_test_PRS_CS <- merge(PRS_test_PRS_CS, data, by=c("IID", "FID"))

  #training
  file_name <- paste(out_PRScs, "/", training_prefix, "_", a, ".profile", sep="")
  data <- read.table(file_name, header=T)
  data <- data[,c(1,2,6)]
  colnames(data) <- c("FID", "IID", b)
  PRS_training_PRS_CS <- merge(PRS_training_PRS_CS, data, by=c("IID", "FID")) 
}

#### Merge files with PCs

PRS_test_PRS_CS_PC <- merge(PRS_test_PRS_CS, PCs_test, by=c("IID", "FID"))
PRS_training_PRS_CS_PC <- merge(PRS_training_PRS_CS, PCs_training, by=c("IID", "FID"))

#### Merge them for PC correction
PRS_test_PRS_CS_PC$Group <- "test"
PRS_training_PRS_CS_PC$Group <- "training"
merged_PRS_CS <- rbind(PRS_test_PRS_CS_PC, PRS_training_PRS_CS_PC)

na_col_PRS_CS <- colSums(is.na(merged_PRS_CS))
all_na_cols_PRS_CS <- which(na_col_PRS_CS == nrow(merged_PRS_CS))
if(length(all_na_cols_PRS_CS) > 0) {
  merged_PRS_CS <- merged_PRS_CS[, -all_na_cols_PRS_CS]
} else {
  # If all_na_cols is empty, just use the original dataframe
  merged_PRS_CS <- merged_PRS_CS
}
na_rows_PRS_CS <- apply(merged_PRS_CS, 1, function(row) any(is.na(row)))
merged_PRS_CS <- merged_PRS_CS[!na_rows_PRS_CS, ]


##############################
#######PC correction

phi_values_res <- paste0(phi_values, "_res")

for (i in 1:length(phi_values)) {
  score <- phi_values[i]
  score_res <- phi_values_res[i]
  merged_PRS_CS <- PC_correction(data=merged_PRS_CS, score=score, score_res=score_res)
}

#### Split again to then do standardization
PRS_test_PRS_CS_PC <- merged_PRS_CS[(merged_PRS_CS$Group=="test"),]
PRS_test_PRS_CS_PC <- PRS_test_PRS_CS_PC[, !colnames(PRS_test_PRS_CS_PC) %in% "Group"]
PRS_training_PRS_CS_PC <- merged_PRS_CS[(merged_PRS_CS$Group=="training"),]
PRS_training_PRS_CS_PC <- PRS_training_PRS_CS_PC[, !colnames(PRS_training_PRS_CS_PC) %in% "Group"]

##############################
#######Standardization

#### Standardize

phi_values_sc <- c(phi_values, phi_values_res)
phi_values_sc_2 <- paste(phi_values_sc, "_sc", sep="")

for (i in 1:length(phi_values_sc)) {
  score <- phi_values_sc[i]
  score_sc <- phi_values_sc_2[i]
  PRS_test_PRS_CS_PC <- standardize(data_ref=PRS_training_PRS_CS_PC, data_target=PRS_test_PRS_CS_PC, score=score, score_sc=score_sc)
  PRS_training_PRS_CS_PC <- standardize(data_ref=PRS_training_PRS_CS_PC, data_target=PRS_training_PRS_CS_PC, score=score, score_sc=score_sc)
}

### Get only scores to write to a file

PRS_test_PRS_CS_scores <- colnames(PRS_test_PRS_CS_PC)[!(colnames(PRS_test_PRS_CS_PC) %in% colnames(PC_only))]
PRS_test_PRS_CS <- PRS_test_PRS_CS_PC[, PRS_test_PRS_CS_scores]
PRS_training_PRS_CS_scores <- colnames(PRS_training_PRS_CS_PC)[!(colnames(PRS_training_PRS_CS_PC) %in% colnames(PC_only))]
PRS_training_PRS_CS <- PRS_training_PRS_CS_PC[, PRS_training_PRS_CS_scores]

write.table(PRS_test_PRS_CS, paste0(out_PRScs, "/", test_prefix, "_scaled_scores"), row.names = F, quote = F, sep="\t")
write.table(PRS_training_PRS_CS, paste0(out_PRScs, "/", training_prefix, "_scaled_scores"), row.names = F, quote = F, sep="\t")


##############################
#######Linear regression

phi_values_res_sc <- paste(phi_values, "_res_sc", sep="")
PRS_test_PRS_CS_pheno <- merge(PRS_test_PRS_CS, pheno, by=c("IID", "FID"))
PRS_test_PRS_CS_pheno$PHENO <- as.numeric(PRS_test_PRS_CS_pheno$PHENO)

#### Perform the regression

if (cross_validation=="FALSE"){
  Regression_results_PRS_CS <- data.frame(matrix(ncol=6, nrow=0))
  
  for (i in 1:length(phi_values_res_sc)) {
    a <- phi_values_res_sc[i]
    thresh <- phi_values[i]  
    tip <- paste("PHENO", "~", a,sep="")
    alpha <- lm(tip, data = PRS_test_PRS_CS_pheno)
    R2 <- summary(alpha)$adj.r.squared
    p <- coef(summary(alpha))[2,4]
    beta <- coef(summary(alpha))[2,1]
    SE <- coef(summary(alpha))[2,2]
    Regression_results_PRS_CS <- rbind(Regression_results_PRS_CS, c("PRS-CS",thresh,beta,SE,p,R2))
    colnames(Regression_results_PRS_CS) <- c("Tool","Parameters", "Beta", "SE", "P_value", "R2")
  }
  
  Regression_results_PRS_CS$R2 <- as.numeric(Regression_results_PRS_CS$R2)
  Regression_results_PRS_CS$P_value <- as.numeric(Regression_results_PRS_CS$P_value)
  Regression_results_PRS_CS$P_value <- ifelse(Regression_results_PRS_CS$P_value==0, 3.4e-314, Regression_results_PRS_CS$P_value)
  
  write.table(Regression_results_PRS_CS, paste0(out_PRScs, "/Regression_results_PRScs"), row.names=F, quote=F, sep="\t")
  
} else if (cross_validation=="TRUE"){
  Regression_results_PRS_CS <- data.frame(matrix(ncol=7, nrow=0))
  train_control <- trainControl(method = "cv", number = 10)
  for (i in 1:length(phi_values_res_sc)) {
    a <- phi_values_res_sc[i]
    thresh <- phi_values[i]  
    tip <- paste("PHENO", "~", a,sep="")
    alpha <- lm(tip, data = PRS_test_PRS_CS_pheno)
    R2 <- summary(alpha)$adj.r.squared
    p <- coef(summary(alpha))[2,4]
    beta <- coef(summary(alpha))[2,1]
    SE <- coef(summary(alpha))[2,2]
    model <- train(as.formula(tip), data=PRS_test_PRS_CS_pheno, method="lm", trControl=train_control, metric="Rsquared")
    R2_CV <- (model$results$Rsquared)
    Regression_results_PRS_CS <- rbind(Regression_results_PRS_CS, c("PRS-CS",thresh,beta,SE,p,R2,R2_CV))
    colnames(Regression_results_PRS_CS) <- c("Tool","Parameters", "Beta", "SE", "P_value", "R2","R2_CV")
  }
  
  Regression_results_PRS_CS$R2 <- as.numeric(Regression_results_PRS_CS$R2)
  Regression_results_PRS_CS$R2_CV <- as.numeric(Regression_results_PRS_CS$R2_CV)
  Regression_results_PRS_CS$P_value <- as.numeric(Regression_results_PRS_CS$P_value)
  Regression_results_PRS_CS$P_value <- ifelse(Regression_results_PRS_CS$P_value==0, 3.4e-314, Regression_results_PRS_CS$P_value)
  
  write.table(Regression_results_PRS_CS, paste0(out_PRScs, "/Regression_results_PRScs"), row.names=F, quote=F, sep="\t")
  
} else {
  "Fill in TRUE or FALSE for cross-validation"
}

#### Make bar plot

if (cross_validation=="FALSE"){
  barplot_gradient(data=Regression_results_PRS_CS, Pt="Parameters", R2="R2", P_value="P_value", all_scores=phi_values, fill_color="cadetblue3", out=paste0(out_PRScs, "/bar_plot_R2_PRScs.svg"), Tool="PRS-CS")
} else if (cross_validation=="TRUE"){
  barplot_gradient(data=Regression_results_PRS_CS, Pt="Parameters", R2="R2_CV", P_value="P_value", all_scores=phi_values, fill_color="cadetblue3", out=paste0(out_PRScs, "/bar_plot_R2_PRScs.svg"), Tool="PRS-CS")
} else {
  "Fill in TRUE or FALSE for cross-validation"
}
 
##############################
#######Select best PRS and plot

if (cross_validation=="FALSE"){
  #### Get best PRS
  Regression_results_valid_PRS_CS <- Regression_results_PRS_CS[(Regression_results_PRS_CS$P_value < 0.05 & Regression_results_PRS_CS$Beta > 0),]
  if (nrow(Regression_results_valid_PRS_CS) == 0) {
    sorted_regression_results_PRS_CS <- data.frame(Tool = "PRS-CS", Parameters = paste0(colnames(PRS_test_PRS_CS[3])), Beta = 0, SE = 0, P_value = 1, R2 = 0.001)
    cat("There are no regression results with p < 0.05 and Beta > 0... creating dummy outcome")
  } else { 
    sorted_regression_results_PRS_CS <- Regression_results_valid_PRS_CS %>% arrange(desc(R2))
  }
  
  best_phi <- sorted_regression_results_PRS_CS[1,2]
  best_phi <- paste0(best_phi, "_res_sc")
  col_select_PRS_CS <- c("FID","IID", best_phi, "PHENO")
  best_PRS_PRS_CS <- select(PRS_test_PRS_CS_pheno, all_of(col_select_PRS_CS))
  best_PRS_PRS_CS$Tool <- "PRS_CS"
  best_PRS_PRS_CS$parameters <- sorted_regression_results_PRS_CS[1,2]
  colnames(best_PRS_PRS_CS) <- c("FID","IID", "score", "PHENO", "Tool", "parameters")
  
  Regression_best_PRS_CS <- sorted_regression_results_PRS_CS[1,]
  Regression_best_per_tool <- rbind(Regression_best_per_tool, Regression_best_PRS_CS)
  
  #### Write to new file
  
  col_select_PRS_CS <- c("FID","IID", "score", "Tool", "parameters")
  selection_best_PRS_PRS_CS <- select(best_PRS_PRS_CS, all_of(col_select_PRS_CS))
  write.table(selection_best_PRS_PRS_CS, paste0(out_PRScs, "/best_PRS_PRScs"), row.names=F, quote=F, sep="\t")
  
} else if (cross_validation=="TRUE"){
  #### Get best PRS
  Regression_results_valid_PRS_CS <- Regression_results_PRS_CS[(Regression_results_PRS_CS$P_value < 0.05 & Regression_results_PRS_CS$Beta > 0),]
  if (nrow(Regression_results_valid_PRS_CS) == 0) {
    sorted_regression_results_PRS_CS <- data.frame(Tool = "PRS-CS", Parameters = paste0(colnames(PRS_test_PRS_CS[3])), Beta = 0, SE = 0, P_value = 1, R2 = 0.001, R2_CV)
    cat("There are no regression results with p < 0.05 and Beta > 0... creating dummy outcome")
  } else { 
    sorted_regression_results_PRS_CS <- Regression_results_valid_PRS_CS %>% arrange(desc(R2_CV))
  }
  
  best_phi <- sorted_regression_results_PRS_CS[1,2]
  best_phi <- paste0(best_phi, "_res_sc")
  col_select_PRS_CS <- c("FID","IID", best_phi, "PHENO")
  best_PRS_PRS_CS <- select(PRS_test_PRS_CS_pheno, all_of(col_select_PRS_CS))
  best_PRS_PRS_CS$Tool <- "PRS_CS"
  best_PRS_PRS_CS$parameters <- sorted_regression_results_PRS_CS[1,2]
  colnames(best_PRS_PRS_CS) <- c("FID","IID", "score", "PHENO", "Tool", "parameters")
  
  Regression_best_PRS_CS <- sorted_regression_results_PRS_CS[1,]
  Regression_best_per_tool <- rbind(Regression_best_per_tool, Regression_best_PRS_CS)
  
  #### Write to new file
  
  col_select_PRS_CS <- c("FID","IID", "score", "Tool", "parameters")
  selection_best_PRS_PRS_CS <- select(best_PRS_PRS_CS, all_of(col_select_PRS_CS))
  write.table(selection_best_PRS_PRS_CS, paste0(out_PRScs, "/best_PRS_PRScs"), row.names=F, quote=F, sep="\t")
  
} else {
  "Fill in TRUE or FALSE for cross-validation"
}

#### Correlation plot

correlation_plot(data=PRS_test_PRS_CS_pheno, x_var=best_phi, y_var="PHENO", point_color="mediumpurple", out=paste0(out_PRScs, "/Correlation_plot_PRScs.svg"))

cat("Done getting best PRS for PRS-CS.... LDpred2....\n")

##################
#####LDpred2######
##################

##############################
#######Load scores

##### Inf
PRS_inf_test <- read.table(paste0(out_LDpred2, "/", "pred_inf_", test_prefix_rsID), header=T)
PRS_inf_training <- read.table(paste0(out_LDpred2, "/", "pred_inf_", training_prefix_rsID), header=T)

##### Grid
PRS_grid_test <- read.table(paste0(out_LDpred2, "/", "pred_grid_", test_prefix_rsID), header=T)
PRS_grid_training <- read.table(paste0(out_LDpred2, "/", "pred_grid_", training_prefix_rsID), header=T)

##### Auto grid
PRS_auto_grid_test <- read.table(paste0(out_LDpred2, "/", "pred_auto_grid_", test_prefix_rsID), header=T)
PRS_auto_grid_training <- read.table(paste0(out_LDpred2, "/", "pred_auto_grid_", training_prefix_rsID), header=T)

##### Auto
PRS_auto_test <- read.table(paste0(out_LDpred2, "/", "pred_auto_", test_prefix_rsID), header=T)
PRS_auto_training <- read.table(paste0(out_LDpred2, "/", "pred_auto_", training_prefix_rsID), header=T)

#### Merge files with PCs

##### Inf
PRS_inf_test_PCs <- merge(PRS_inf_test, PCs_test, by=c("IID", "FID"))
PRS_inf_training_PCs <- merge(PRS_inf_training, PCs_training, by=c("IID", "FID"))

##### Grid
PRS_grid_test_PCs <- merge(PRS_grid_test, PCs_test, by=c("IID", "FID"))
PRS_grid_training_PCs <- merge(PRS_grid_training, PCs_training, by=c("IID", "FID"))

##### Auto grid
PRS_auto_grid_test_PCs <- merge(PRS_auto_grid_test, PCs_test, by=c("IID", "FID"))
PRS_auto_grid_training_PCs <- merge(PRS_auto_grid_training, PCs_training, by=c("IID", "FID"))

##### Auto
PRS_auto_test_PCs <- merge(PRS_auto_test, PCs_test, by=c("IID", "FID"))
PRS_auto_training_PCs <- merge(PRS_auto_training, PCs_training, by=c("IID", "FID"))

##############################
#######PC correction

#### Merge test and training for PC correction

##### Inf
PRS_inf_test_PCs$Group <- "test"
PRS_inf_training_PCs$Group <- "training"
merge_inf <- rbind(PRS_inf_test_PCs, PRS_inf_training_PCs)
merge_inf <- merge_inf[!is.na(merge_inf$pred_inf), ]

##### Grid
PRS_grid_test_PCs$Group <- "test"
PRS_grid_training_PCs$Group <- "training"
merge_grid <- rbind(PRS_grid_test_PCs, PRS_grid_training_PCs)
na_col_grid <- colSums(is.na(merge_grid))
all_na_cols_grid <- which(na_col_grid == nrow(merge_grid))
if(length(all_na_cols_grid) > 0) {
  merge_grid <- merge_grid[, -all_na_cols_grid]
} else {
  # If all_na_cols is empty, just use the original dataframe
  merge_grid <- merge_grid
}
na_rows_grid <- apply(merge_grid, 1, function(row) any(is.na(row)))
merge_grid <- merge_grid[!na_rows_grid, ]

##### Auto grid
PRS_auto_grid_test_PCs$Group <- "test"
PRS_auto_grid_training_PCs$Group <- "training"
merge_auto_grid <- rbind(PRS_auto_grid_test_PCs, PRS_auto_grid_training_PCs)
na_col_auto_grid <- colSums(is.na(merge_auto_grid))
all_na_cols_auto_grid <- which(na_col_auto_grid == nrow(merge_auto_grid))
if(length(all_na_cols_auto_grid) > 0) {
  merge_auto_grid <- merge_auto_grid[, -all_na_cols_auto_grid]
} else {
  # If all_na_cols is empty, just use the original dataframe
  merge_auto_grid <- merge_auto_grid
}
na_rows_auto_grid <- apply(merge_auto_grid, 1, function(row) any(is.na(row)))
merge_auto_grid <- merge_auto_grid[!na_rows_auto_grid, ]

##### Auto
PRS_auto_test_PCs$Group <- "test"
PRS_auto_training_PCs$Group <- "training"
merge_auto <- rbind(PRS_auto_test_PCs, PRS_auto_training_PCs)
merge_auto <- merge_auto[!is.na(merge_auto$pred_auto),]

####PC correction

##### Inf
merge_inf <- PC_correction(data=merge_inf, score="pred_inf", score_res="pred_inf_res")

##### Grid
scores_grid <- select(merge_grid, -IID, -FID, -Group, -colnames(PC_only))
scores_grid_2 <- colnames(scores_grid)
scores_grid_res <- paste0(scores_grid_2, "_res")
for (i in 1:length(scores_grid_2)){
  score <- scores_grid_2[i]
  score_res <- scores_grid_res[i]
  merge_grid <- PC_correction(data=merge_grid, score=score, score_res=score_res)
}

##### Auto grid
scores_auto_grid <- select(merge_auto_grid, -IID, -FID, -Group, -colnames(PC_only))
scores_auto_grid_2 <- colnames(scores_auto_grid)
scores_auto_grid_res <- paste0(scores_auto_grid_2, "_res")
for (i in 1:length(scores_auto_grid_2)){
  score <- scores_auto_grid_2[i]
  score_res <- scores_auto_grid_res[i]
  merge_auto_grid <- PC_correction(data=merge_auto_grid, score=score, score_res=score_res)
}

##### Auto
merge_auto <- PC_correction(data=merge_auto, score="pred_auto", score_res = "pred_auto_res")


#### Split again to then do standardization

##### Inf
PRS_inf_test_PCs <- merge_inf[(merge_inf$Group=="test"),]
PRS_inf_test_PCs <- PRS_inf_test_PCs[, !colnames(PRS_inf_test_PCs) %in% "Group"]
PRS_inf_training_PCs <- merge_inf[(merge_inf$Group=="training"),]
PRS_inf_training_PCs <- PRS_inf_training_PCs[, !colnames(PRS_inf_training_PCs) %in% "Group"]

##### Grid
PRS_grid_test_PCs <- merge_grid[(merge_grid$Group=="test"),]
PRS_grid_test_PCs <- PRS_grid_test_PCs[, !colnames(PRS_grid_test_PCs) %in% "Group"]
PRS_grid_training_PCs <- merge_grid[(merge_grid$Group=="training"),]
PRS_grid_training_PCs <- PRS_grid_training_PCs[, !colnames(PRS_grid_training_PCs) %in% "Group"]

##### Auto Grid
PRS_auto_grid_test_PCs <- merge_auto_grid[(merge_auto_grid$Group=="test"),]
PRS_auto_grid_test_PCs <- PRS_auto_grid_test_PCs[, !colnames(PRS_auto_grid_test_PCs) %in% "Group"]
PRS_auto_grid_training_PCs <- merge_auto_grid[(merge_auto_grid$Group=="training"),]
PRS_auto_grid_training_PCs <- PRS_auto_grid_training_PCs[, !colnames(PRS_auto_grid_training_PCs) %in% "Group"]

##### Auto
PRS_auto_test_PCs <- merge_auto[(merge_auto$Group=="test"),]
PRS_auto_test_PCs <- PRS_auto_test_PCs[, !colnames(PRS_auto_test_PCs) %in% "Group"]
PRS_auto_training_PCs <- merge_auto[(merge_auto$Group=="training"),]
PRS_auto_training_PCs <- PRS_auto_training_PCs[, !colnames(PRS_auto_training_PCs) %in% "Group"]

##############################
#######Standardization

#### Standardize

##### Inf
sc_inf <- c("pred_inf", "pred_inf_res")
sc2_inf <- paste0(sc_inf, "_sc")
for (i in 1:length(sc_inf)){
  score <- sc_inf[i]
  score_sc <- sc2_inf[i]
  PRS_inf_training_PCs <- standardize(data_ref=PRS_inf_training_PCs, data_target = PRS_inf_training_PCs, score=score, score_sc=score_sc)
  PRS_inf_test_PCs <- standardize(data_ref=PRS_inf_training_PCs, data_target = PRS_inf_test_PCs, score=score, score_sc=score_sc)
}

##### Grid
sc_grid <- c(scores_grid_2, scores_grid_res)
sc2_grid <- paste0(sc_grid, "_sc")
for (i in 1:length(sc_grid)){
  score <- sc_grid[i]
  score_sc <- sc2_grid[i]
  PRS_grid_training_PCs <- standardize(data_ref=PRS_grid_training_PCs, data_target = PRS_grid_training_PCs, score=score, score_sc=score_sc)
  PRS_grid_test_PCs <- standardize(data_ref=PRS_grid_training_PCs, data_target = PRS_grid_test_PCs, score=score, score_sc=score_sc)
}

##### Auto grid
sc_auto_grid <- c(scores_auto_grid_2, scores_auto_grid_res)
sc2_auto_grid <- paste0(sc_auto_grid, "_sc")
for (i in 1:length(sc_auto_grid)){
  score <- sc_auto_grid[i]
  score_sc <- sc2_auto_grid[i]
  PRS_auto_grid_training_PCs <- standardize(data_ref=PRS_auto_grid_training_PCs, data_target = PRS_auto_grid_training_PCs, score=score, score_sc=score_sc)
  PRS_auto_grid_test_PCs <- standardize(data_ref=PRS_auto_grid_training_PCs, data_target = PRS_auto_grid_test_PCs, score=score, score_sc=score_sc)
}

##### Auto
sc_auto <- c("pred_auto", "pred_auto_res")
sc2_auto <- paste0(sc_auto, "_sc")
for (i in 1:length(sc_auto)){
  score <- sc_auto[i]
  score_sc <- sc2_auto[i]
  PRS_auto_training_PCs <- standardize(data_ref=PRS_auto_training_PCs, data_target = PRS_auto_training_PCs, score=score, score_sc=score_sc)
  PRS_auto_test_PCs <- standardize(data_ref=PRS_auto_training_PCs, data_target = PRS_auto_test_PCs, score=score, score_sc=score_sc)
}

#### Get only scores to write to a file

##### Inf
PRS_inf_test_scores <- colnames(PRS_inf_test_PCs)[!(colnames(PRS_inf_test_PCs) %in% colnames(PC_only))]
PRS_inf_test <- PRS_inf_test_PCs[, PRS_inf_test_scores]
PRS_inf_training_scores <- colnames(PRS_inf_training_PCs)[!(colnames(PRS_inf_training_PCs) %in% colnames(PC_only))]
PRS_inf_training <- PRS_inf_training_PCs[, PRS_inf_training_scores]
write.table(PRS_inf_test, paste0(out_LDpred2, "/", test_prefix_rsID, "_scaled_scores_inf"), row.names = F, quote = F, sep="\t")
write.table(PRS_inf_training, paste0(out_LDpred2, "/", training_prefix_rsID, "_scaled_scores_inf"), row.names = F, quote = F, sep="\t")


##### Grid
PRS_grid_test_scores <- colnames(PRS_grid_test_PCs)[!(colnames(PRS_grid_test_PCs) %in% colnames(PC_only))]
PRS_grid_test <- PRS_grid_test_PCs[, PRS_grid_test_scores]
PRS_grid_training_scores <- colnames(PRS_grid_training_PCs)[!(colnames(PRS_grid_training_PCs) %in% colnames(PC_only))]
PRS_grid_training <- PRS_grid_training_PCs[, PRS_grid_training_scores]
write.table(PRS_grid_test, paste0(out_LDpred2, "/", test_prefix_rsID, "_scaled_scores_grid"), row.names = F, quote = F, sep="\t")
write.table(PRS_grid_training, paste0(out_LDpred2, "/", training_prefix_rsID, "_scaled_scores_grid"), row.names = F, quote = F, sep="\t")

##### Auto grid
PRS_auto_grid_test_scores <- colnames(PRS_auto_grid_test_PCs)[!(colnames(PRS_auto_grid_test_PCs) %in% colnames(PC_only))]
PRS_auto_grid_test <- PRS_auto_grid_test_PCs[, PRS_auto_grid_test_scores]
PRS_auto_grid_training_scores <- colnames(PRS_auto_grid_training_PCs)[!(colnames(PRS_auto_grid_training_PCs) %in% colnames(PC_only))]
PRS_auto_grid_training <- PRS_auto_grid_training_PCs[, PRS_auto_grid_training_scores]
write.table(PRS_auto_grid_test, paste0(out_LDpred2, "/", test_prefix_rsID, "_scaled_scores_auto_grid"), row.names = F, quote = F, sep="\t")
write.table(PRS_auto_grid_training, paste0(out_LDpred2, "/", training_prefix_rsID, "_scaled_scores_auto_grid"), row.names = F, quote = F, sep="\t")


##### Auto
PRS_auto_test_scores <- colnames(PRS_auto_test_PCs)[!(colnames(PRS_auto_test_PCs) %in% colnames(PC_only))]
PRS_auto_test <- PRS_auto_test_PCs[, PRS_auto_test_scores]
PRS_auto_training_scores <- colnames(PRS_auto_training_PCs)[!(colnames(PRS_auto_training_PCs) %in% colnames(PC_only))]
PRS_auto_training <- PRS_auto_training_PCs[, PRS_auto_training_scores]
write.table(PRS_auto_test, paste0(out_LDpred2, "/", test_prefix_rsID, "_scaled_scores_auto"), row.names = F, quote = F, sep="\t")
write.table(PRS_auto_training, paste0(out_LDpred2, "/", training_prefix_rsID, "_scaled_scores_auto"), row.names = F, quote = F, sep="\t")


##############################
#######Linear regression

#### Prepare files for regression

##### Inf
scores_regr_inf <- colnames(PRS_inf_test)[grep("_res_sc$", colnames(PRS_inf_test))]
columns_inf <- c("IID", "FID", scores_regr_inf)
PRS_inf_for_regression <- PRS_inf_test[, columns_inf]
PRS_inf_for_regression <- merge(PRS_inf_for_regression, pheno, by=c("IID", "FID"))
colnames(PRS_inf_for_regression) <- c("IID", "FID",scores_regr_inf, "PHENO")
PRS_inf_for_regression$PHENO <- as.numeric(PRS_inf_for_regression$PHENO)

##### Grid
scores_regr_grid <- colnames(PRS_grid_test)[grep("_res_sc$", colnames(PRS_grid_test))]
columns_grid <- c("IID", "FID", scores_regr_grid)
PRS_grid_for_regression <- PRS_grid_test[, columns_grid]
PRS_grid_for_regression <- merge(PRS_grid_for_regression, pheno, by=c("IID", "FID"))
colnames(PRS_grid_for_regression) <- c("IID", "FID",scores_regr_grid, "PHENO")
PRS_grid_for_regression$PHENO <- as.numeric(PRS_grid_for_regression$PHENO)

##### Auto grid
scores_regr_auto_grid <- colnames(PRS_auto_grid_test)[grep("_res_sc$", colnames(PRS_auto_grid_test))]
columns_auto_grid <- c("IID", "FID", scores_regr_auto_grid)
PRS_auto_grid_for_regression <- PRS_auto_grid_test[, columns_auto_grid]
PRS_auto_grid_for_regression <- merge(PRS_auto_grid_for_regression, pheno, by=c("IID", "FID"))
colnames(PRS_auto_grid_for_regression) <- c("IID", "FID",scores_regr_auto_grid, "PHENO")
PRS_auto_grid_for_regression$PHENO <- as.numeric(PRS_auto_grid_for_regression$PHENO)

##### Auto
scores_regr_auto <- colnames(PRS_auto_test)[grep("_res_sc$", colnames(PRS_auto_test))]
columns_auto <- c("IID", "FID", scores_regr_auto)
PRS_auto_for_regression <- PRS_auto_test[, columns_auto]
PRS_auto_for_regression <- merge(PRS_auto_for_regression, pheno, by=c("IID", "FID"))
colnames(PRS_auto_for_regression) <- c("IID", "FID",scores_regr_auto, "PHENO")
PRS_auto_for_regression$PHENO <- as.numeric(PRS_auto_for_regression$PHENO)

#### Perform linear regression

##### Inf

if (cross_validation=="FALSE"){
  Regression_results_inf <- data.frame(matrix(ncol=7, nrow=0))
  for (i in 1:length(scores_regr_inf)) {
    b <- scores_regr_inf[i]
    thresh <- paste(b)  
    tip <- paste("PHENO", "~", b,sep="")
    alpha <- lm(tip, data = PRS_inf_for_regression)
    R2 <- summary(alpha)$adj.r.squared
    p <- coef(summary(alpha))[2,4]
    beta <- coef(summary(alpha))[2,1]
    SE <- coef(summary(alpha))[2,2]
    Regression_results_inf <- rbind(Regression_results_inf, c("inf", thresh,beta,SE,p,R2, "LDpred2"))
    colnames(Regression_results_inf) <- c("Model", "Parameters", "Beta", "SE", "P_value", "R2", "Tool")
  }
  Regression_results_inf$R2 <- as.numeric(Regression_results_inf$R2)
  Regression_results_inf$P_value <- as.numeric(Regression_results_inf$P_value)
  Regression_results_inf$P_value <- ifelse(Regression_results_inf$P_value==0, 3.4e-314, Regression_results_inf$P_value)
  
  write.table(Regression_results_inf, paste0(out_LDpred2, "/Regression_results_LDpred2_inf"), row.names=F, quote=F, sep="\t")
  
} else if (cross_validation=="TRUE"){
  Regression_results_inf <- data.frame(matrix(ncol=8, nrow=0))
  train_control <- trainControl(method = "cv", number = 10)
  for (i in 1:length(scores_regr_inf)) {
    b <- scores_regr_inf[i]
    thresh <- paste(b)  
    tip <- paste("PHENO", "~", b,sep="")
    alpha <- lm(tip, data = PRS_inf_for_regression)
    R2 <- summary(alpha)$adj.r.squared
    p <- coef(summary(alpha))[2,4]
    beta <- coef(summary(alpha))[2,1]
    SE <- coef(summary(alpha))[2,2]
    model <- train(as.formula(tip), data=PRS_inf_for_regression, method="lm", trControl=train_control, metric="Rsquared")
    R2_CV <- (model$results$Rsquared)
    Regression_results_inf <- rbind(Regression_results_inf, c("inf", thresh,beta,SE,p,R2,R2_CV, "LDpred2"))
    colnames(Regression_results_inf) <- c("Model", "Parameters", "Beta", "SE", "P_value", "R2","R2_CV", "Tool")
  }
  Regression_results_inf$R2 <- as.numeric(Regression_results_inf$R2)
  Regression_results_inf$R2_CV <- as.numeric(Regression_results_inf$R2_CV)
  Regression_results_inf$P_value <- as.numeric(Regression_results_inf$P_value)
  Regression_results_inf$P_value <- ifelse(Regression_results_inf$P_value==0, 3.4e-314, Regression_results_inf$P_value)
  
  write.table(Regression_results_inf, paste0(out_LDpred2, "/Regression_results_LDpred2_inf"), row.names=F, quote=F, sep="\t")
  
} else {
  "Fill in TRUE or FALSE for cross-validation"
}


##### Grid

if (cross_validation=="FALSE"){
  Regression_results_grid <- data.frame(matrix(ncol=7, nrow=0))
  for (i in 1:length(scores_regr_grid)) {
    b <- scores_regr_grid[i]
    thresh <- paste(b)  
    tip <- paste("PHENO", "~", b,sep="")
    alpha <- lm(tip, data = PRS_grid_for_regression)
    R2 <- summary(alpha)$adj.r.squared
    p <- coef(summary(alpha))[2,4]
    beta <- coef(summary(alpha))[2,1]
    SE <- coef(summary(alpha))[2,2]
    Regression_results_grid <- rbind(Regression_results_grid, c("grid", thresh,beta,SE,p,R2, "LDpred2"))
    colnames(Regression_results_grid) <- c("Model", "Parameters", "Beta", "SE", "P_value", "R2", "Tool")
  }
  
  Regression_results_grid$R2 <- as.numeric(Regression_results_grid$R2)
  Regression_results_grid$P_value <- as.numeric(Regression_results_grid$P_value)
  Regression_results_grid$P_value <- ifelse(Regression_results_grid$P_value==0, 3.4e-314, Regression_results_grid$P_value)
  
  write.table(Regression_results_grid, paste0(out_LDpred2, "/Regression_results_LDpred2_grid"), row.names=F, quote=F, sep="\t")
  
} else if (cross_validation=="TRUE"){
  Regression_results_grid <- data.frame(matrix(ncol=8, nrow=0))
  train_control <- trainControl(method = "cv", number = 10)
  for (i in 1:length(scores_regr_grid)) {
    b <- scores_regr_grid[i]
    thresh <- paste(b)  
    tip <- paste("PHENO", "~", b,sep="")
    alpha <- lm(tip, data = PRS_grid_for_regression)
    R2 <- summary(alpha)$adj.r.squared
    p <- coef(summary(alpha))[2,4]
    beta <- coef(summary(alpha))[2,1]
    SE <- coef(summary(alpha))[2,2]
    model <- train(as.formula(tip), data=PRS_grid_for_regression, method="lm", trControl=train_control, metric="Rsquared")
    R2_CV <- (model$results$Rsquared)
    Regression_results_grid <- rbind(Regression_results_grid, c("grid", thresh,beta,SE,p,R2, R2_CV,"LDpred2"))
    colnames(Regression_results_grid) <- c("Model", "Parameters", "Beta", "SE", "P_value", "R2", "R2_CV","Tool")
  }
  
  Regression_results_grid$R2 <- as.numeric(Regression_results_grid$R2)
  Regression_results_grid$R2_CV <- as.numeric(Regression_results_grid$R2_CV)
  Regression_results_grid$P_value <- as.numeric(Regression_results_grid$P_value)
  Regression_results_grid$P_value <- ifelse(Regression_results_grid$P_value==0, 3.4e-314, Regression_results_grid$P_value)
  
  write.table(Regression_results_grid, paste0(out_LDpred2, "/Regression_results_LDpred2_grid"), row.names=F, quote=F, sep="\t")
  
} else {
  "Fill in TRUE or FALSE for cross-validation"
}


##### Auto grid

if (cross_validation=="FALSE"){
  Regression_results_auto_grid <- data.frame(matrix(ncol=7, nrow=0))
  for (i in 1:length(scores_regr_auto_grid)) {
    b <- scores_regr_auto_grid[i]
    thresh <- paste(b)  
    tip <- paste("PHENO", "~", b,sep="")
    alpha <- lm(tip, data = PRS_auto_grid_for_regression)
    R2 <- summary(alpha)$adj.r.squared
    p <- coef(summary(alpha))[2,4]
    beta <- coef(summary(alpha))[2,1]
    SE <- coef(summary(alpha))[2,2]
    Regression_results_auto_grid <- rbind(Regression_results_auto_grid, c("auto_grid", thresh,beta,SE,p,R2, "LDpred2"))
    colnames(Regression_results_auto_grid) <- c("Model", "Parameters", "Beta", "SE", "P_value", "R2", "Tool")
  }
  
  Regression_results_auto_grid$R2 <- as.numeric(Regression_results_auto_grid$R2)
  Regression_results_auto_grid$P_value <- as.numeric(Regression_results_auto_grid$P_value)
  Regression_results_auto_grid$P_value <- ifelse(Regression_results_auto_grid$P_value==0, 3.4e-314, Regression_results_auto_grid$P_value)
  
  write.table(Regression_results_auto_grid, paste0(out_LDpred2, "/Regression_results_LDpred2_auto_grid"), row.names=F, quote=F, sep="\t")
  
} else if (cross_validation=="TRUE"){
  Regression_results_auto_grid <- data.frame(matrix(ncol=8, nrow=0))
  train_control <- trainControl(method = "cv", number = 10)
  for (i in 1:length(scores_regr_auto_grid)) {
    b <- scores_regr_auto_grid[i]
    thresh <- paste(b)  
    tip <- paste("PHENO", "~", b,sep="")
    alpha <- lm(tip, data = PRS_auto_grid_for_regression)
    R2 <- summary(alpha)$adj.r.squared
    p <- coef(summary(alpha))[2,4]
    beta <- coef(summary(alpha))[2,1]
    SE <- coef(summary(alpha))[2,2]
    model <- train(as.formula(tip), data=PRS_auto_grid_for_regression, method="lm", trControl=train_control, metric="Rsquared")
    R2_CV <- (model$results$Rsquared)
    Regression_results_auto_grid <- rbind(Regression_results_auto_grid, c("auto_grid", thresh,beta,SE,p,R2,R2_CV, "LDpred2"))
    colnames(Regression_results_auto_grid) <- c("Model", "Parameters", "Beta", "SE", "P_value", "R2", "R2_CV","Tool")
  }
  
  Regression_results_auto_grid$R2 <- as.numeric(Regression_results_auto_grid$R2)
  Regression_results_auto_grid$R2_CV <- as.numeric(Regression_results_auto_grid$R2_CV)
  Regression_results_auto_grid$P_value <- as.numeric(Regression_results_auto_grid$P_value)
  Regression_results_auto_grid$P_value <- ifelse(Regression_results_auto_grid$P_value==0, 3.4e-314, Regression_results_auto_grid$P_value)
  
  write.table(Regression_results_auto_grid, paste0(out_LDpred2, "/Regression_results_LDpred2_auto_grid"), row.names=F, quote=F, sep="\t")
  
} else {
  "Fill in TRUE or FALSE for cross-validation"
}

##### Auto

if (cross_validation=="FALSE"){
  Regression_results_auto <- data.frame(matrix(ncol=7, nrow=0))
  for (i in 1:length(scores_regr_auto)) {
    b <- scores_regr_auto[i]
    thresh <- paste(b)  
    tip <- paste("PHENO", "~", b,sep="")
    alpha <- lm(tip, data = PRS_auto_for_regression)
    R2 <- summary(alpha)$adj.r.squared
    p <- coef(summary(alpha))[2,4]
    beta <- coef(summary(alpha))[2,1]
    SE <- coef(summary(alpha))[2,2]
    Regression_results_auto <- rbind(Regression_results_auto, c("auto", thresh,beta,SE,p,R2, "LDpred2"))
    colnames(Regression_results_auto) <- c("Model", "Parameters", "Beta", "SE", "P_value", "R2", "Tool")
  }
  
  Regression_results_auto$R2 <- as.numeric(Regression_results_auto$R2)
  Regression_results_auto$P_value <- as.numeric(Regression_results_auto$P_value)
  Regression_results_auto$P_value <- ifelse(Regression_results_auto$P_value==0, 3.4e-314, Regression_results_auto$P_value)
  
  write.table(Regression_results_auto, paste0(out_LDpred2, "/Regression_results_LDpred2_auto"), row.names=F, quote=F, sep="\t")
  
} else if (cross_validation=="TRUE"){
  Regression_results_auto <- data.frame(matrix(ncol=8, nrow=0))
  train_control <- trainControl(method = "cv", number = 10)
  for (i in 1:length(scores_regr_auto)) {
    b <- scores_regr_auto[i]
    thresh <- paste(b)  
    tip <- paste("PHENO", "~", b,sep="")
    alpha <- lm(tip, data = PRS_auto_for_regression)
    R2 <- summary(alpha)$adj.r.squared
    p <- coef(summary(alpha))[2,4]
    beta <- coef(summary(alpha))[2,1]
    SE <- coef(summary(alpha))[2,2]
    model <- train(as.formula(tip), data=PRS_auto_for_regression, method="lm", trControl=train_control, metric="Rsquared")
    R2_CV <- (model$results$Rsquared)
    Regression_results_auto <- rbind(Regression_results_auto, c("auto", thresh,beta,SE,p,R2, R2_CV,"LDpred2"))
    colnames(Regression_results_auto) <- c("Model", "Parameters", "Beta", "SE", "P_value", "R2", "R2_CV","Tool")
  }
  
  Regression_results_auto$R2 <- as.numeric(Regression_results_auto$R2)
  Regression_results_auto$R2_CV <- as.numeric(Regression_results_auto$R2_CV)
  Regression_results_auto$P_value <- as.numeric(Regression_results_auto$P_value)
  Regression_results_auto$P_value <- ifelse(Regression_results_auto$P_value==0, 3.4e-314, Regression_results_auto$P_value)
  
  write.table(Regression_results_auto, paste0(out_LDpred2, "/Regression_results_LDpred2_auto"), row.names=F, quote=F, sep="\t")
  
} else {
  "Fill in TRUE or FALSE for cross-validation"
}


#### Make bar plot

##### Inf

if (cross_validation=="FALSE"){
  bar_plot <- ggplot(data=Regression_results_inf, aes(x=Parameters, y=R2)) +
    geom_bar(stat = "identity", position = "dodge", fill="cadetblue3") +
    scale_x_discrete(limits=scores_regr_inf, guide = guide_axis(angle = 90)) +
    xlab("") + ylab("R²") +theme_minimal() + ggtitle("LDpred2 - inf model")
  ggsave(paste0(out_LDpred2, "/bar_plot_R2_LDpred2_inf.svg"), bar_plot, width=7, height=7)
} else if (cross_validation=="TRUE"){
  bar_plot <- ggplot(data=Regression_results_inf, aes(x=Parameters, y=R2_CV)) +
    geom_bar(stat = "identity", position = "dodge", fill="cadetblue3") +
    scale_x_discrete(limits=scores_regr_inf, guide = guide_axis(angle = 90)) +
    xlab("") + ylab("R²") +theme_minimal() + ggtitle("LDpred2 - inf model")
  ggsave(paste0(out_LDpred2, "/bar_plot_R2_LDpred2_inf.svg"), bar_plot, width=7, height=7)
} else {
  "Fill in TRUE or FALSE for cross-validation"
}

##### Grid

if (cross_validation=="FALSE"){
  barplot_gradient(data=Regression_results_grid, Pt="Parameters", R2="R2", P_value="P_value", all_scores=scores_regr_grid, fill_color="cadetblue3", out=paste0(out_LDpred2, "/bar_plot_R2_LDpred2_grid.svg"), Tool="LDpred2 - grid model")
} else if (cross_validation=="TRUE"){
  barplot_gradient(data=Regression_results_grid, Pt="Parameters", R2="R2_CV", P_value="P_value", all_scores=scores_regr_grid, fill_color="cadetblue3", out=paste0(out_LDpred2, "/bar_plot_R2_LDpred2_grid.svg"), Tool="LDpred2 - grid model")
} else {
  "Fill in TRUE or FALSE for cross-validation"
}

##### Auto grid

if (cross_validation=="FALSE"){
  barplot_gradient(data=Regression_results_auto_grid, Pt="Parameters", R2="R2", P_value="P_value", all_scores=scores_regr_auto_grid, fill_color="cadetblue3", out=paste0(out_LDpred2, "/bar_plot_R2_LDpred2_auto_grid.svg"), Tool="LDpred2 - auto grid model")
} else if (cross_validation=="TRUE"){
  barplot_gradient(data=Regression_results_auto_grid, Pt="Parameters", R2="R2_CV", P_value="P_value", all_scores=scores_regr_auto_grid, fill_color="cadetblue3", out=paste0(out_LDpred2, "/bar_plot_R2_LDpred2_auto_grid.svg"), Tool="LDpred2 - auto grid model")
} else {
  "Fill in TRUE or FALSE for cross-validation"
}

##### Auto

if (cross_validation=="FALSE"){
  bar_plot <- ggplot(data=Regression_results_auto, aes(x=Parameters, y=R2)) +
    geom_bar(stat = "identity", position = "dodge", fill="cadetblue3") +
    scale_x_discrete(limits=scores_regr_auto, guide = guide_axis(angle = 90)) +
    xlab("") + ylab("Nagelkerke R2") +theme_minimal() + ggtitle("LDpred2 - auto model")
  ggsave(paste0(out_LDpred2, "/bar_plot_R2_LDpred2_auto.svg"), bar_plot, width=7, height=7)
} else if (cross_validation=="TRUE"){
  bar_plot <- ggplot(data=Regression_results_auto, aes(x=Parameters, y=R2_CV)) +
    geom_bar(stat = "identity", position = "dodge", fill="cadetblue3") +
    scale_x_discrete(limits=scores_regr_auto, guide = guide_axis(angle = 90)) +
    xlab("") + ylab("Nagelkerke R2") +theme_minimal() + ggtitle("LDpred2 - auto model")
  ggsave(paste0(out_LDpred2, "/bar_plot_R2_LDpred2_auto.svg"), bar_plot, width=7, height=7)
} else {
  "Fill in TRUE or FALSE for cross-validation"
}

##############################
#######Select best PRS and plot

#### Get best PRS
Regression_results_LDpred2 <- rbind(Regression_results_inf, rbind(Regression_results_grid, rbind(Regression_results_auto_grid, Regression_results_auto)))
Regression_results_valid_LDpred2 <- Regression_results_LDpred2[(Regression_results_LDpred2$P_value < 0.05 & Regression_results_LDpred2$Beta > 0),]

if (cross_validation=="FALSE"){
  if (nrow(Regression_results_valid_LDpred2) == 0) {
    sorted_regression_results_LDpred2 <- data.frame(Model = "grid", Parameters = paste0(colnames(PRS_grid_training[3]), "_res_sc"), Beta = 0, SE = 0, P_value = 1, R2 = 0.001, Tool="LDpred2")
    cat("There are no regression results with p < 0.05 and Beta > 0... creating dummy outcome")
  } else { 
    sorted_regression_results_LDpred2 <- Regression_results_valid_LDpred2 %>% arrange(desc(R2))
  }
  
  best_model <- sorted_regression_results_LDpred2[1,1]
  best_parameters <- sorted_regression_results_LDpred2[1,2]
  
  best_model_LDpred2 <- switch(
    best_model,
    "inf" = PRS_inf_for_regression,
    "grid" = PRS_grid_for_regression,
    "auto_grid" = PRS_auto_grid_for_regression,
    "auto" = PRS_auto_for_regression,
    stop("Invalid best_model value")
  )
  
  col_select_LDpred2 <- c("FID","IID", best_parameters, "PHENO")
  best_PRS_LDpred2 <- select(best_model_LDpred2, all_of(col_select_LDpred2))
  best_PRS_LDpred2$Tool <- paste0("LDpred2_", best_model)
  best_PRS_LDpred2$parameters <- best_parameters
  colnames(best_PRS_LDpred2) <- c("FID","IID", "score", "PHENO", "Tool", "parameters")
  
  #### Write to new file
  
  col_select_LDpred2 <- c("FID","IID", "score", "Tool", "parameters")
  selection_best_PRS_LDpred2 <- select(best_PRS_LDpred2, all_of(col_select_LDpred2))
  write.table(selection_best_PRS_LDpred2, paste0(out_LDpred2, "/best_PRS_LDpred2"), row.names=F, quote=F, sep="\t")
  
  Regression_best_LDpred2 <- sorted_regression_results_LDpred2[1,]
  Regression_best_LDpred2$Parameters <- paste0(Regression_best_LDpred2$Model, "_", Regression_best_LDpred2$Parameters)
  Regression_best_LDpred2 <- Regression_best_LDpred2[,c(2:7)]
  Regression_best_per_tool <- rbind(Regression_best_per_tool, Regression_best_LDpred2)
} else if (cross_validation=="TRUE"){
  if (nrow(Regression_results_valid_LDpred2) == 0) {
    sorted_regression_results_LDpred2 <- data.frame(Model = "grid", Parameters = paste0(colnames(PRS_grid_training[3]), "_res_sc"), Beta = 0, SE = 0, P_value = 1, R2 = 0.001, R2_CV = 0.001,Tool="LDpred2")
    cat("There are no regression results with p < 0.05 and Beta > 0... creating dummy outcome")
  } else { 
    sorted_regression_results_LDpred2 <- Regression_results_valid_LDpred2 %>% arrange(desc(R2_CV))
  }
  
  best_model <- sorted_regression_results_LDpred2[1,1]
  best_parameters <- sorted_regression_results_LDpred2[1,2]
  
  best_model_LDpred2 <- switch(
    best_model,
    "inf" = PRS_inf_for_regression,
    "grid" = PRS_grid_for_regression,
    "auto_grid" = PRS_auto_grid_for_regression,
    "auto" = PRS_auto_for_regression,
    stop("Invalid best_model value")
  )
  
  col_select_LDpred2 <- c("FID","IID", best_parameters, "PHENO")
  best_PRS_LDpred2 <- select(best_model_LDpred2, all_of(col_select_LDpred2))
  best_PRS_LDpred2$Tool <- paste0("LDpred2_", best_model)
  best_PRS_LDpred2$parameters <- best_parameters
  colnames(best_PRS_LDpred2) <- c("FID","IID", "score", "PHENO", "Tool", "parameters")
  
  #### Write to new file
  
  col_select_LDpred2 <- c("FID","IID", "score", "Tool", "parameters")
  selection_best_PRS_LDpred2 <- select(best_PRS_LDpred2, all_of(col_select_LDpred2))
  write.table(selection_best_PRS_LDpred2, paste0(out_LDpred2, "/best_PRS_LDpred2"), row.names=F, quote=F, sep="\t")
  
  Regression_best_LDpred2 <- sorted_regression_results_LDpred2[1,]
  Regression_best_LDpred2$Parameters <- paste0(Regression_best_LDpred2$Model, "_", Regression_best_LDpred2$Parameters)
  Regression_best_LDpred2 <- Regression_best_LDpred2[,c(2:8)]
  Regression_best_per_tool <- rbind(Regression_best_per_tool, Regression_best_LDpred2)
} else {
  "Fill in TRUE or FALSE for cross-validation"
}

#### Correlation plot

correlation_plot(data=best_PRS_LDpred2, x_var="score", y_var="PHENO", point_color="mediumpurple", out=paste0(out_LDpred2, "/Correlation_plot_LDpred2.svg"))

cat("Done getting best PRS for LDpred2.... Lassosum2....\n")

##################
####Lassosum-2####
##################

##############################
#######Load scores

PRS_lasso2_test <- read.table(paste0(out_lasso2, "/", "pred_lasso2_", test_prefix_rsID), header=T)
PRS_lasso2_training <- read.table(paste0(out_lasso2, "/", "pred_lasso2_", training_prefix_rsID), header=T)

#### Merge files with PCs

PRS_lasso2_test_PCs <- merge(PRS_lasso2_test, PCs_test, by=c("IID", "FID"))
PRS_lasso2_training_PCs <- merge(PRS_lasso2_training, PCs_training, by=c("IID", "FID"))

##############################
#######PC correction

#### Merge test and training for PC correction

PRS_lasso2_test_PCs$Group <- "test"
PRS_lasso2_training_PCs$Group <- "training"
merge_lasso2 <- rbind(PRS_lasso2_test_PCs, PRS_lasso2_training_PCs)
na_col_lasso2 <- colSums(is.na(merge_lasso2))
all_na_cols_lasso2 <- which(na_col_lasso2 == nrow(merge_lasso2))
if(length(all_na_cols_lasso2) > 0) {
  merge_lasso2 <- merge_lasso2[, -all_na_cols_lasso2]
} else {
  # If all_na_cols is empty, just use the original dataframe
  merge_lasso2 <- merge_lasso2
}
na_rows_lasso2 <- apply(merge_lasso2, 1, function(row) any(is.na(row)))
merge_lasso2 <- na.omit(merge_lasso2)

####PC correction

scores_lasso2 <- select(merge_lasso2, -IID, -FID, -Group, -colnames(PC_only))
scores_lasso2_2 <- colnames(scores_lasso2)
scores_lasso2_res <- paste0(scores_lasso2_2, "_res")
for (i in 1:length(scores_lasso2_2)){
  score <- scores_lasso2_2[i]
  score_res <- scores_lasso2_res[i]
  merge_lasso2 <- PC_correction(data=merge_lasso2, score=score, score_res=score_res)
}

#### Split again to then do standardization

PRS_lasso2_test_PCs <- merge_lasso2[(merge_lasso2$Group=="test"),]
PRS_lasso2_test_PCs <- PRS_lasso2_test_PCs[, !colnames(PRS_lasso2_test_PCs) %in% "Group"]
PRS_lasso2_training_PCs <- merge_lasso2[(merge_lasso2$Group=="training"),]
PRS_lasso2_training_PCs <- PRS_lasso2_training_PCs[, !colnames(PRS_lasso2_training_PCs) %in% "Group"]

##############################
#######Standardization

#### Standardize

sc_lasso2 <- c(scores_lasso2_2, scores_lasso2_res)
sc2_lasso2 <- paste0(sc_lasso2, "_sc")
for (i in 1:length(sc_lasso2)){
  score <- sc_lasso2[i]
  score_sc <- sc2_lasso2[i]
  PRS_lasso2_training_PCs <- standardize(data_ref=PRS_lasso2_training_PCs, data_target = PRS_lasso2_training_PCs, score=score, score_sc=score_sc)
  PRS_lasso2_test_PCs <- standardize(data_ref=PRS_lasso2_training_PCs, data_target = PRS_lasso2_test_PCs, score=score, score_sc=score_sc)
}

#### Get only scores to write to a file

PRS_lasso2_test_scores <- colnames(PRS_lasso2_test_PCs)[!(colnames(PRS_lasso2_test_PCs) %in% colnames(PC_only))]
PRS_lasso2_test <- PRS_lasso2_test_PCs[, PRS_lasso2_test_scores]
PRS_lasso2_training_scores <- colnames(PRS_lasso2_training_PCs)[!(colnames(PRS_lasso2_training_PCs) %in% colnames(PC_only))]
PRS_lasso2_training <- PRS_lasso2_training_PCs[, PRS_lasso2_training_scores]
write.table(PRS_lasso2_test, paste0(out_lasso2, "/", test_prefix_rsID, "_scaled_scores_lasso2"), row.names = F, quote = F, sep="\t")
write.table(PRS_lasso2_training, paste0(out_lasso2, "/", training_prefix_rsID, "_scaled_scores_lasso2"), row.names = F, quote = F, sep="\t")

##############################
#######Linear regression

#### Prepare files for regression
scores_regr_lasso2 <- colnames(PRS_lasso2_test)[grep("_res_sc$", colnames(PRS_lasso2_test))]
columns_lasso2 <- c("IID", "FID", scores_regr_lasso2)
PRS_lasso2_for_regression <- PRS_lasso2_test[, columns_lasso2]
PRS_lasso2_for_regression <- merge(PRS_lasso2_for_regression, pheno, by=c("IID", "FID"))
colnames(PRS_lasso2_for_regression) <- c("IID", "FID",scores_regr_lasso2, "PHENO")
PRS_lasso2_for_regression$PHENO <- as.numeric(PRS_lasso2_for_regression$PHENO)

#### Perform linear regression

if (cross_validation=="FALSE"){
  Regression_results_lasso2 <- data.frame(matrix(ncol=6, nrow=0))
  for (i in 1:length(scores_regr_lasso2)) {
    b <- scores_regr_lasso2[i]
    thresh <- paste(b)  
    tip <- paste("PHENO", "~", b,sep="")
    alpha <- lm(tip, data = PRS_lasso2_for_regression)
    R2 <- summary(alpha)$adj.r.squared
    p <- coef(summary(alpha))[2,4]
    beta <- coef(summary(alpha))[2,1]
    SE <- coef(summary(alpha))[2,2]
    Regression_results_lasso2 <- rbind(Regression_results_lasso2, c("lassosum2", thresh,beta,SE,p,R2))
    colnames(Regression_results_lasso2) <- c("Tool", "Parameters", "Beta", "SE", "P_value", "R2")
  }
  
  Regression_results_lasso2$R2 <- as.numeric(Regression_results_lasso2$R2)
  Regression_results_lasso2$P_value <- as.numeric(Regression_results_lasso2$P_value)
  Regression_results_lasso2$P_value <- ifelse(Regression_results_lasso2$P_value==0, 3.4e-314, Regression_results_lasso2$P_value)
  
  write.table(Regression_results_lasso2, paste0(out_lasso2, "/Regression_results_lassosum2"), row.names=F, quote=F, sep="\t")
} else if (cross_validation=="TRUE"){
  Regression_results_lasso2 <- data.frame(matrix(ncol=7, nrow=0))
  train_control <- trainControl(method = "cv", number = 10)
  for (i in 1:length(scores_regr_lasso2)) {
    b <- scores_regr_lasso2[i]
    thresh <- paste(b)  
    tip <- paste("PHENO", "~", b,sep="")
    alpha <- lm(tip, data = PRS_lasso2_for_regression)
    R2 <- summary(alpha)$adj.r.squared
    p <- coef(summary(alpha))[2,4]
    beta <- coef(summary(alpha))[2,1]
    SE <- coef(summary(alpha))[2,2]
    model <- train(as.formula(tip), data=PRS_lasso2_for_regression, method="lm", trControl=train_control, metric="Rsquared")
    R2_CV <- (model$results$Rsquared)
    Regression_results_lasso2 <- rbind(Regression_results_lasso2, c("lassosum2", thresh,beta,SE,p,R2, R2_CV))
    colnames(Regression_results_lasso2) <- c("Tool", "Parameters", "Beta", "SE", "P_value", "R2", "R2_CV")
  }
  
  Regression_results_lasso2$R2 <- as.numeric(Regression_results_lasso2$R2)
  Regression_results_lasso2$R2_CV <- as.numeric(Regression_results_lasso2$R2_CV)
  Regression_results_lasso2$P_value <- as.numeric(Regression_results_lasso2$P_value)
  Regression_results_lasso2$P_value <- ifelse(Regression_results_lasso2$P_value==0, 3.4e-314, Regression_results_lasso2$P_value)
  
  write.table(Regression_results_lasso2, paste0(out_lasso2, "/Regression_results_lassosum2"), row.names=F, quote=F, sep="\t")
} else {
  "Fill in TRUE or FALSE for cross-validation"
}

#### Make bar plot

if (cross_validation=="FALSE"){
  barplot_gradient(data=Regression_results_lasso2, Pt="Parameters", R2="R2", P_value="P_value", all_scores=scores_regr_lasso2, fill_color="cadetblue3", out=paste0(out_LDpred2, "/bar_plot_R2_lassosum2.svg"), Tool="Lassosum2")
} else if (cross_validation=="TRUE"){
  barplot_gradient(data=Regression_results_lasso2, Pt="Parameters", R2="R2_CV", P_value="P_value", all_scores=scores_regr_lasso2, fill_color="cadetblue3", out=paste0(out_LDpred2, "/bar_plot_R2_lassosum2.svg"), Tool="Lassosum2")
} else {
  "Fill in TRUE or FALSE for cross-validation"
}

##############################
#######Select best PRS and plot

#### Get best PRS
Regression_results_valid_lasso2 <- Regression_results_lasso2[(Regression_results_lasso2$P_value < 0.05 & Regression_results_lasso2$Beta > 0),]

if (cross_validation=="FALSE"){
  if (nrow(Regression_results_valid_lasso2) == 0) {
    sorted_regression_results_lasso2 <- data.frame(Tool = "lassosum2", Parameters=paste0(colnames(PRS_lasso2_test)[3], "_res_sc"), Beta = 0, SE = 0, P_value = 1, R2 = 0.001)
    cat("There are no regression results with p < 0.05 and Beta > 0... creating dummy outcome")
  } else { 
    sorted_regression_results_lasso2 <- Regression_results_valid_lasso2 %>% arrange(desc(R2))
  }
  
  best_model <- sorted_regression_results_lasso2[1,1]
  best_parameters <- sorted_regression_results_lasso2[1,2]
  col_select_lasso2 <- c("FID","IID", best_parameters, "PHENO")
  best_PRS_lasso2 <- select(PRS_lasso2_for_regression, all_of(col_select_lasso2))
  best_PRS_lasso2$Tool <- best_model
  best_PRS_lasso2$parameters <- best_parameters
  colnames(best_PRS_lasso2) <- c("FID","IID", "score", "PHENO", "Tool", "parameters")
  
  #### Write to new file
  
  col_select_lasso2 <- c("FID","IID", "score", "Tool", "parameters")
  selection_best_PRS_lasso2 <- select(best_PRS_lasso2, all_of(col_select_lasso2))
  write.table(selection_best_PRS_lasso2, paste0(out_lasso2, "/best_PRS_lassosum2"), row.names=F, quote=F, sep="\t")
  
  Regression_best_lasso2 <- sorted_regression_results_lasso2[1,]
  Regression_best_per_tool <- rbind(Regression_best_per_tool, Regression_best_lasso2)
  write.table(Regression_best_per_tool, paste0(out_comparison, "/Regression_results_best_per_tool"), quote = F, row.names = F ,sep="\t")
  
} else if (cross_validation=="TRUE"){
  if (nrow(Regression_results_valid_lasso2) == 0) {
    sorted_regression_results_lasso2 <- data.frame(Tool = "lassosum2", Parameters=paste0(colnames(PRS_lasso2_test)[3], "_res_sc"), Beta = 0, SE = 0, P_value = 1, R2 = 0.001, R2_CV=0.001)
    cat("There are no regression results with p < 0.05 and Beta > 0... creating dummy outcome")
  } else { 
    sorted_regression_results_lasso2 <- Regression_results_valid_lasso2 %>% arrange(desc(R2_CV))
  }
  
  best_model <- sorted_regression_results_lasso2[1,1]
  best_parameters <- sorted_regression_results_lasso2[1,2]
  col_select_lasso2 <- c("FID","IID", best_parameters, "PHENO")
  best_PRS_lasso2 <- select(PRS_lasso2_for_regression, all_of(col_select_lasso2))
  best_PRS_lasso2$Tool <- best_model
  best_PRS_lasso2$parameters <- best_parameters
  colnames(best_PRS_lasso2) <- c("FID","IID", "score", "PHENO", "Tool", "parameters")
  
  #### Write to new file
  
  col_select_lasso2 <- c("FID","IID", "score", "Tool", "parameters")
  selection_best_PRS_lasso2 <- select(best_PRS_lasso2, all_of(col_select_lasso2))
  write.table(selection_best_PRS_lasso2, paste0(out_lasso2, "/best_PRS_lassosum2"), row.names=F, quote=F, sep="\t")
  
  Regression_best_lasso2 <- sorted_regression_results_lasso2[1,]
  Regression_best_per_tool <- rbind(Regression_best_per_tool, Regression_best_lasso2)
  write.table(Regression_best_per_tool, paste0(out_comparison, "/Regression_results_best_per_tool"), quote = F, row.names = F ,sep="\t")
  
} else {
  "Fill in TRUE or FALSE for cross-validation"
}

#### Correlation plot

correlation_plot(data=best_PRS_lasso2, x_var="score", y_var="PHENO", point_color="mediumpurple", out=paste0(out_lasso2, "/Correlation_plot_lassosum2.svg"))

cat("Done getting best PRS for lassosum2...\n")
cat("Done getting best PRS per tool!\nGet best PRS overall...\n")

##################
####Comparison####
##################

#### Compare R2

if (cross_validation=="FALSE"){
  Regression_best_per_tool$Tool <- factor(Regression_best_per_tool$Tool, levels = c("PRSice", "PRS-CS", "LDpred2","lassosum", "lassosum2"))
  Comparison_bar_plot <- ggplot(data=Regression_best_per_tool, aes(x=Tool, y=R2, fill=Tool)) +
    geom_bar(stat="identity", position="dodge") +
    scale_fill_manual(values=c("#CC79A7", "#56B4E9","cyan3", "#E69F00","coral"), name="Tool") +
    theme_minimal()+
    theme(legend.position = "none") +
    theme(axis.text.x = element_text(size=12, face = "bold"),
          axis.text.y = element_text(size=12, face = "bold"),
          axis.title.y = element_text(size=15,face="bold"),
          axis.title.x = element_text(size=15,face="bold")) +
    xlab("") + ylab("R²\n") + ggtitle("Comparison R² per tool")
  ggsave(paste0(out_comparison, "/Bar_plot_comparison_R2_per_tool.svg"), width=7, height=5)
} else if (cross_validation=="TRUE"){
  Regression_best_per_tool$Tool <- factor(Regression_best_per_tool$Tool, levels = c("PRSice", "PRS-CS", "LDpred2","lassosum", "lassosum2"))
  Comparison_bar_plot <- ggplot(data=Regression_best_per_tool, aes(x=Tool, y=R2_CV, fill=Tool)) +
    geom_bar(stat="identity", position="dodge") +
    scale_fill_manual(values=c("#CC79A7", "#56B4E9","cyan3", "#E69F00","coral"), name="Tool") +
    theme_minimal()+
    theme(legend.position = "none") +
    theme(axis.text.x = element_text(size=12, face = "bold"),
          axis.text.y = element_text(size=12, face = "bold"),
          axis.title.y = element_text(size=15,face="bold"),
          axis.title.x = element_text(size=15,face="bold")) +
    xlab("") + ylab("R²\n") + ggtitle("Comparison R² per tool")
  ggsave(paste0(out_comparison, "/Bar_plot_comparison_R2_per_tool.svg"), width=7, height=5)
  
} else {
  "Fill in TRUE or FALSE for cross-validation"
}

cat("Done comparing tools!\nPipeline completed!\n")
