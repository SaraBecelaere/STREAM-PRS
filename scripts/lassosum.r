#!/usr/bin/env Rscript

#Load required packages
rm(list = ls())
args = commandArgs(trailingOnly=TRUE)

library(dplyr)
library(data.table)
library(methods)
library(magrittr)
library(devtools)
library(parallel)
library(remotes)

list.of.packages <- c("lassosum")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[, "Package"])]
if(length(new.packages)){
  remotes::install_version("Matrix", version="1.6.5", repos='http://cran.us.r-project.org')
  install_github("tshmak/lassosum")
} 
library(lassosum)

set.seed(17)

test_prefix = args[1]
training_prefix = args[13]
bfile_test <- args[2]
sum.stat <- args[3]
bfile_training <- args[4]
ld.blocks <- args[5]
chromosome_col <- args[6]
pos_col <- args[7]
effect_allele <- args[8]
noneffect_allele <- args[9]
N_col <- args[10]
lasso_thresholding_values <- args[11]
lambda <- args[14]
out_lasso <- args[12]

#Adapt some parameters
lasso_thresholding_values <- strsplit(lasso_thresholding_values, ",")[[1]]
lasso_thresholding_values <- as.numeric(lasso_thresholding_values)
lambda <- eval(parse(text=lambda))

#Read in the summary statistics
ss <- fread(sum.stat)

#Change N that are too small to minimum allowed N
ss$N <- pull(ss, N_col)
ss$N <- as.numeric(ss$N)
minimum=max(ss$N, 30)/10
ss$N[ss$N < minimum] <- minimum

#Transform the P-values into correlation
cor <- p2cor(p = ss$P_num, n = ss$N, sign = ss$beta_num)

#Lassosum for training data

##obtain fam file
fam <- fread(paste0(bfile_training, ".fam"))
fam[,ID:=do.call(paste, c(.SD, sep=":")),.SDcols=c(1:2)]

##Run the lassosum pipeline
out <- lassosum.pipeline(
  cor = cor,
  chr = pull(ss, chromosome_col),
  pos = pull(ss, pos_col),
  A1 = pull(ss, effect_allele),
  A2 = pull(ss, noneffect_allele),
  ref.bfile = bfile_training,
  test.bfile = bfile_training,
  LDblocks = ld.blocks,
  lambda=lambda,
  s=lasso_thresholding_values
)

##Lassosum for test data
out2 <- subset(out, s=out$s, lambda=out$lambda)
dummy_pheno <- rnorm(nrow.bfile(bfile_test))
v2 <- validate(out2, test.bfile=bfile_test, plot=FALSE, pheno=dummy_pheno)

###extract scores and betas for every lambda per s and write to a file

names <- c(sprintf("L_%f", lambda))
fam_test <- fread(paste0(bfile_test, ".fam"))
fam_test[,ID:=do.call(paste, c(.SD, sep=":")),.SDcols=c(1:2)]

for (i in 1:length(lasso_thresholding_values)){
  s_value <- lasso_thresholding_values[i]
  data_frame <- as.data.frame(out$pgs[[i]])
  colnames(data_frame) <- names
  data_frame$FID <- fam$V1
  data_frame$IID <- fam$V2
  write.table(data_frame, quote = FALSE, row.names = FALSE, col.names= TRUE, paste0(out_lasso, "/", training_prefix,"_scores_s_", s_value, ".txt"))
  
  data_frame_test <- as.data.frame(v2$pgs[[i]])
  colnames(data_frame_test) <- names
  data_frame_test$FID <- fam_test$V1
  data_frame_test$IID <- fam_test$V2
  write.table(data_frame_test, quote = FALSE, row.names = FALSE, col.names= TRUE, paste0(out_lasso, "/", test_prefix,"_scores_s_", s_value, ".txt"))
  
  data_frame_2 <- as.data.frame(out$beta[[i]])
  colnames(data_frame_2) <- names
  write.table(data_frame_2, quote = FALSE, row.names = FALSE, col.names = TRUE , paste0(out_lasso, "/", training_prefix,"_betas_", s_value, ".txt"))
}

###extract SNPs info
results_snp <- out$sumstats
write.table(results_snp, quote = FALSE, row.names = FALSE, col.names = FALSE, paste0(out_lasso, "/", training_prefix,"_snps.txt"))

cat("Lassosum completed!\n")


