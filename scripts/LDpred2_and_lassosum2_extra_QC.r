#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

library(remotes)
list.of.packages <- c("bigsnpr")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos='http://cran.us.r-project.org')
library(bigsnpr)
library(bigreadr)
library(data.table)
library(magrittr)
options(bigstatsr.check.parallel.blas = FALSE)
options(default.nproc.blas = NULL)
library(dplyr)
library(ggplot2)

set.seed(564)

ss <- args[1]
training_prefix <- args[2]
training_file <- args[3]
test_prefix <- args[4]
test_file <- args[5]
SNP_col <- args[6]
CHR_col <- args[7]
POS_col <- args[8]
effect_A <- args[9]
other_A <- args[10]
N_col <- args[11]
beta_se_col <- args[12]
allele_freq_col <- args[13]
out_LDpred2 <- args[14]
out_lasso2 <- args[15]
ldref_hm3_plus <- args[16]
NCORES <- args[17]
NCORES <- as.numeric(NCORES)

#Grid model parameters
N_h2 <- args[18]
N_h2 <- eval(parse(text=N_h2))
N_p <- args[19]
N_p <- eval(parse(text=N_p))

#Auto grid model parameters
initial_p <- args[20]
initial_p <- eval(parse(text=initial_p))

#Lassosum2 parameters
delta <- args[21]
delta <- eval(parse(text=delta))
nlambda <- args[22]
nlambda <- as.numeric(nlambda)
min_ratio <- args[23]
min_ratio <- as.numeric(min_ratio)

map_path <- args[24]

cat("Starting LDpred2...\n")

map <- readRDS(paste0(map_path))

sumstats <- bigreadr::fread2(ss)
sumstats$chr <- as.numeric(pull(sumstats, CHR_col))
sumstats$pos <- pull(sumstats, POS_col)
sumstats$a1 <- pull(sumstats, effect_A)
sumstats$a0 <- pull(sumstats, other_A)
sumstats$freq <- pull(sumstats, allele_freq_col)
sumstats$n_eff <- pull(sumstats, N_col)
sumstats$beta <- sumstats$beta_num
sumstats$beta_se <- pull(sumstats, beta_se_col)
sumstats$rsid <- pull(sumstats, SNP_col)

# Matching and QC

info_snp <- as_tibble(snp_match(sumstats, map, return_flip_and_rev = TRUE, join_by_pos=FALSE))

info_snp2 <- info_snp %>%
  filter(n_eff > (0.8 * max(n_eff)), !`_FLIP_`, !`_REV_`) %>%
  mutate(sd_af = sqrt(2 * freq * (1 - freq)),
         sd_ss = 1 / sqrt(n_eff * beta_se^2 + beta^2),
         sd_ss = sd_ss * print(sqrt(0.5) / quantile(sd_ss, 0.99)))

info_snp2$is_bad <- with(info_snp2,
                         sd_ss < (0.7 * sd_af) | sd_ss > (sd_af + 0.1) |
                           sd_ss < 0.05 | sd_af < 0.05)

df_beta <- filter(info_snp2, !is_bad)
#df_beta <- info_snp
df_beta2 <- data.frame(df_beta)

cat("SNPs matched! Building correlation matrix...\n")

# Build correlation matrix

tmp <- tempfile(tmpdir = paste0(out_LDpred2, "/tmp-data"))
for (chr in 1:22) {
  
  cat(chr, ".. ", sep = "")
  
  ## indices in 'df_beta'
  ind.chr <- which(df_beta$chr == chr)
  ## indices in 'map_ldref'
  ind.chr2 <- df_beta$`_NUM_ID_`[ind.chr]
  ## indices in 'corr_chr'
  ind.chr3 <- match(ind.chr2, which(map$chr == chr))
  
  corr_chr <- readRDS(paste0(ldref_hm3_plus, "/LD_with_blocks_chr", chr, ".rds"))[ind.chr3, ind.chr3]
  
  if (chr == 1) {
    corr <- as_SFBM(corr_chr, tmp, compact = TRUE)
  } else {
    corr$add_columns(corr_chr, nrow(corr))
  }
}
corr

cat("Correlation matrix built! Estimating heritability...\n")

# Estimate heritability

(ldsc <- snp_ldsc2(corr, df_beta, blocks = 200, intercept = NULL, ncores = NCORES))

ldsc_h2_est <- ldsc[["h2"]]

if (ldsc_h2_est < 0) {
  ldsc_h2_est <- 0.001
  cat("Heritability estimate is negative, changed value to 0.001\n")
}


cat("Heritability estimated! Loading in training and test data...\n")

# Training data

if (!file.exists(paste0(training_file, ".rds"))) {
  # If the file doesn't exist, create a new file
  cat("training .rds file doesn't exist. Creating it now...")
  snp_readBed(paste0(training_file, ".bed"))
} else {
  # If the file already exists, do nothing or print a message
  cat("training .rds file already exists. No new file created.\n")
}

obj.bigSNP_training <- snp_attach(paste0(training_file, ".rds"))

G_training   <- obj.bigSNP_training$genotypes
G_training_imp <- snp_fastImputeSimple(G_training, method = "mean2", ncores = NCORES)
G_training <- G_training_imp

sample_ids_training <- obj.bigSNP_training$fam$sample.ID
map_bigSNP_training <- dplyr::transmute(obj.bigSNP_training$map,
                                        chr = as.integer(chromosome), pos = physical.pos,
                                        a0 = allele1, a1 = allele2, rsid=marker.ID)

map_pgs_training <- df_beta[1:4]; map_pgs_training$beta <- 1
map_pgs2_training <- snp_match(map_pgs_training, map_bigSNP_training, join_by_pos=FALSE)

# Test data

if (!file.exists(paste0(test_file, ".rds"))) {
  # If the file doesn't exist, create a new file
  cat("test .rds file doesn't exist. Creating it now...")
  snp_readBed(paste0(test_file, ".bed"))
} else {
  # If the file already exists, do nothing or print a message
  cat("test .rds file already exists. No new file created.\n")
}

obj.bigSNP_test <- snp_attach(paste0(test_file, ".rds"))

G_test   <- obj.bigSNP_test$genotypes
G_test_imp <- snp_fastImputeSimple(G_test, method = "mean2", ncores = NCORES)
G_test <- G_test_imp

sample_ids_test <- obj.bigSNP_test$fam$sample.ID
map_bigSNP_test <- dplyr::transmute(obj.bigSNP_test$map,
                                    chr = as.integer(chromosome), pos = physical.pos,
                                    a0 = allele1, a1 = allele2, rsid=marker.ID)

map_pgs_test <- snp_match(map_pgs2_training, map_bigSNP_test, join_by_pos=FALSE)

cat("Data loaded and SNPs matched! Running infinitesimal model...\n")

#### INF model

beta_inf <- snp_ldpred2_inf(corr, df_beta, h2 = ldsc_h2_est)
pred_inf_training <- big_prodVec(G_training, beta_inf[map_pgs2_training[["_NUM_ID_.ss"]]], ind.col = map_pgs2_training[["_NUM_ID_"]])
pred_inf_test <- big_prodVec(G_test, beta_inf[map_pgs_test[["_NUM_ID_.ss"]]], ind.col = map_pgs_test[["_NUM_ID_"]])

beta_inf2 <- as.data.frame(beta_inf)
beta_inf2$rsid <- df_beta2$rsid

pred_inf_training2 <- as.data.frame(pred_inf_training)
colnames(pred_inf_training2) <- "pred_inf"
pred_inf_training2$FID <- sample_ids_training
pred_inf_training2$IID <- sample_ids_training

pred_inf_test2 <- as.data.frame(pred_inf_test)
colnames(pred_inf_test2) <- "pred_inf"
pred_inf_test2$FID <- sample_ids_test
pred_inf_test2$IID <- sample_ids_test


write.table(beta_inf2, paste0(out_LDpred2, "/beta_inf"), row.names = F, quote = F, sep="\t")
write.table(pred_inf_training2, paste0(out_LDpred2, "/pred_inf_", training_prefix), row.names = F, quote = F, sep="\t")
write.table(pred_inf_test2, paste0(out_LDpred2, "/pred_inf_", test_prefix), row.names = F, quote = F, sep="\t")


cat("Infinitesimal model completed! Running grid model...\n")

#### GRID model
h2_seq <- round(ldsc_h2_est * N_h2, 4)
p_seq <- signif(N_p, 2)
params <- expand.grid(p = p_seq, h2 = h2_seq, sparse = c(FALSE, TRUE))

beta_grid <- snp_ldpred2_grid(corr, df_beta, params)
pred_grid_training <- big_prodMat(G_training, beta_grid[map_pgs2_training[["_NUM_ID_.ss"]], ], ind.col = map_pgs2_training[["_NUM_ID_"]])
pred_grid_test <- big_prodMat(G_test, beta_grid[map_pgs_test[["_NUM_ID_.ss"]], ], ind.col = map_pgs_test[["_NUM_ID_"]])

params$colname <- paste0("p_",params$p,"_h2_",params$h2,"_sparse_",params$sparse)

beta_grid_withcolnames <- as.data.frame(beta_grid)
colnames_betagrid <- as.character(params[,4])
colnames(beta_grid_withcolnames) <- colnames_betagrid
beta_grid_withcolnames$rsid <- df_beta2$rsid

pred_grid_withcolnames_training <- as.data.frame(pred_grid_training)
colnames(pred_grid_withcolnames_training) <- colnames_betagrid
pred_grid_withcolnames_training$FID <- sample_ids_training
pred_grid_withcolnames_training$IID <- sample_ids_training

pred_grid_withcolnames_test <- as.data.frame(pred_grid_test)
colnames(pred_grid_withcolnames_test) <- colnames_betagrid
pred_grid_withcolnames_test$FID <- sample_ids_test
pred_grid_withcolnames_test$IID <- sample_ids_test

write.table(beta_grid_withcolnames, paste0(out_LDpred2, "/beta_grid"), row.names = F, quote = F, sep="\t")
write.table(pred_grid_withcolnames_training, paste0(out_LDpred2, "/pred_grid_", training_prefix), row.names = F, quote = F, sep="\t")
write.table(pred_grid_withcolnames_test, paste0(out_LDpred2, "/pred_grid_", test_prefix), row.names = F, quote = F, sep="\t")


cat("Grid model completed! Running auto grid model...\n")

#### AUTO GRID model

coef_shrink <- 0.95
multi_auto <- snp_ldpred2_auto(corr, df_beta, h2_init = ldsc_h2_est,
                               vec_p_init = initial_p,
                               allow_jump_sign = FALSE, shrink_corr = coef_shrink,
                               ncores = NCORES)
auto <- multi_auto[[1]]

range <- sapply(multi_auto, function(auto) diff(range(auto$corr_est)))
keep <- which(range > (0.95 * quantile(range, 0.95, na.rm = TRUE)))

beta_auto_grid <- sapply(multi_auto[keep], function(auto) auto$beta_est)
pred_auto_grid_training <- big_prodMat(G_training, beta_auto_grid[map_pgs2_training[["_NUM_ID_.ss"]], ], ind.col = map_pgs2_training[["_NUM_ID_"]])
pred_auto_grid_test <- big_prodMat(G_test, beta_auto_grid[map_pgs_test[["_NUM_ID_.ss"]], ], ind.col = map_pgs_test[["_NUM_ID_"]])

colnames_auto <- initial_p[keep]
colnames_auto <- as.data.frame(colnames_auto)
colnames_auto$colname <- paste0("p_init_",colnames_auto$colnames_auto)

beta_auto_grid_withcolnames <- as.data.frame(beta_auto_grid)
colnames_beta_auto_grid <- as.character(colnames_auto[,2])
colnames(beta_auto_grid_withcolnames) <- colnames_beta_auto_grid
beta_auto_grid_withcolnames$rsid <- df_beta2$rsid

pred_auto_grid_withcolnames_training <- as.data.frame(pred_auto_grid_training)
colnames(pred_auto_grid_withcolnames_training) <- colnames_beta_auto_grid
pred_auto_grid_withcolnames_training$FID <- sample_ids_training
pred_auto_grid_withcolnames_training$IID <- sample_ids_training

pred_auto_grid_withcolnames_test <- as.data.frame(pred_auto_grid_test)
colnames(pred_auto_grid_withcolnames_test) <- colnames_beta_auto_grid
pred_auto_grid_withcolnames_test$FID <- sample_ids_test
pred_auto_grid_withcolnames_test$IID <- sample_ids_test

write.table(beta_auto_grid_withcolnames, paste0(out_LDpred2, "/beta_auto_grid"), row.names = F, quote = F, sep="\t")
write.table(pred_auto_grid_withcolnames_training, paste0(out_LDpred2,"/pred_auto_grid_", training_prefix), row.names = F, quote = F, sep="\t")
write.table(pred_auto_grid_withcolnames_test, paste0(out_LDpred2, "/pred_auto_grid_", test_prefix), row.names = F, quote = F, sep="\t")

cat("Auto grid model completed! Running auto model...\n")

#### AUTO model

beta_auto <- rowMeans(sapply(multi_auto[keep], function(auto) auto$beta_est))
pred_auto_training <- big_prodVec(G_training, beta_auto[map_pgs2_training[["_NUM_ID_.ss"]]], ind.col = map_pgs2_training[["_NUM_ID_"]])
pred_auto_test <- big_prodVec(G_test, beta_auto[map_pgs_test[["_NUM_ID_.ss"]]], ind.col = map_pgs_test[["_NUM_ID_"]])

beta_auto2 <- as.data.frame(beta_auto)
beta_auto2$rsid <- df_beta2$rsid

pred_auto2_training <- as.data.frame(pred_auto_training)
colnames(pred_auto2_training) <- "pred_auto"
pred_auto2_training$FID <- sample_ids_training
pred_auto2_training$IID <- sample_ids_training

pred_auto2_test <- as.data.frame(pred_auto_test)
colnames(pred_auto2_test) <- "pred_auto"
pred_auto2_test$FID <- sample_ids_test
pred_auto2_test$IID <- sample_ids_test

write.table(beta_auto2, paste0(out_LDpred2, "/beta_auto"), row.names = F, quote = F, sep="\t")
write.table(pred_auto2_training, paste0(out_LDpred2, "/pred_auto_", training_prefix), row.names = F, quote = F, sep="\t")
write.table(pred_auto2_test, paste0(out_LDpred2, "/pred_auto_", test_prefix), row.names = F, quote = F, sep="\t")

cat("Auto model completed!\n LDpred2 completed!\n Running lassosum2...\n")

#### LASSOSUM 2

beta_lassosum2 <- snp_lassosum2(corr, df_beta, delta = delta, nlambda = nlambda, lambda.min.ratio = min_ratio)
pred_lassosum2_training <- big_prodMat(G_training, beta_lassosum2[map_pgs2_training[["_NUM_ID_.ss"]], ], ind.col = map_pgs2_training[["_NUM_ID_"]])
pred_lassosum2_test <- big_prodMat(G_test, beta_lassosum2[map_pgs_test[["_NUM_ID_.ss"]], ], ind.col = map_pgs_test[["_NUM_ID_"]])

params_lasso <- attr(beta_lassosum2, "grid_param")
params_lasso$colname <- paste0("delta_",params_lasso$delta,"_lambda_",params_lasso$lambda)

beta_lasso_withcolnames <- as.data.frame(beta_lassosum2)
colnames_betalasso <- as.character(params_lasso[,6])
colnames(beta_lasso_withcolnames) <- colnames_betalasso
beta_lasso_withcolnames$rsid <- df_beta2$rsid

pred_lasso_withcolnames_training <- as.data.frame(pred_lassosum2_training)
colnames(pred_lasso_withcolnames_training) <- colnames_betalasso
pred_lasso_withcolnames_training$FID <- sample_ids_training
pred_lasso_withcolnames_training$IID <- sample_ids_training

pred_lasso_withcolnames_test <- as.data.frame(pred_lassosum2_test)
colnames(pred_lasso_withcolnames_test) <- colnames_betalasso
pred_lasso_withcolnames_test$FID <- sample_ids_test
pred_lasso_withcolnames_test$IID <- sample_ids_test

write.table(beta_lasso_withcolnames, paste0(out_lasso2, "/beta_lasso2"), row.names = F, quote = F, sep="\t")
write.table(pred_lasso_withcolnames_training, paste0(out_lasso2, "/pred_lasso2_", training_prefix), row.names = F, quote = F, sep="\t")
write.table(pred_lasso_withcolnames_test, paste0(out_lasso2, "/pred_lasso2_", test_prefix), row.names = F, quote = F, sep="\t")

cat("Lassosum2 completed!\n")


