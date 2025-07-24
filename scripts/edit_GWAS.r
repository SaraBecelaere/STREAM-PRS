#!/usr/bin/env Rscript

args = commandArgs(trailingOnly=TRUE)

library(dplyr)

cat("Start editing GWAS....\n")

sum.stat <- args[1]
sum.stat.2 <- args[2]
p_value_col <- args[3]
beta_col <- args[4]
chr_pos_ID_col <- args[5]
rsID_col <- args[6]
effect_allele <- args[7] 
noneffect_allele <- args[8]
position_col <- args[9]

ss <- read.table(sum.stat, header=T)

ss$P_num <- as.numeric(pull(ss, p_value_col))
ss$P_num[ss$P_num <= 5e-308] <- 5e-308
ss$beta_num <- as.numeric(pull(ss, beta_col))

ss <- ss[(nchar(as.character(pull(ss, effect_allele)))==1 & nchar(as.character(pull(ss, noneffect_allele)))==1),]

ss <- ss[!((as.character(pull(ss, effect_allele))=="A" & as.character(pull(ss, noneffect_allele))=="T") | 
               (as.character(pull(ss, effect_allele))=="T" & as.character(pull(ss, noneffect_allele))=="A") | 
               (as.character(pull(ss, effect_allele))=="C" & as.character(pull(ss, noneffect_allele))=="G") |
               (as.character(pull(ss, effect_allele))=="G" & as.character(pull(ss, noneffect_allele))=="C")),]

duplicates <- ss[duplicated(pull(ss, chr_pos_ID_col)), 1]
ss <- ss[!(pull(ss, chr_pos_ID_col)) %in% duplicates, ]

ss <- ss %>% mutate(across(all_of(position_col), ~ format(., scientific = FALSE)))

write.table(ss, paste0(sum.stat, "_edited_by_pipeline"), row.names=F, quote=F, sep="\t")

cat("Edited GWAS for PRSice and lassosum....\n")

#Prepare PRS-CS and LDpred2 GWAS file

ss_2 <- read.table(sum.stat.2, header=T)

ss_2$P_num <- as.numeric(pull(ss_2, p_value_col))
ss_2$P_num[ss_2$P_num <= 5e-308] <- 5e-308
ss_2$beta_num <- as.numeric(pull(ss_2, beta_col))

ss_2 <- ss_2[(nchar(as.character(pull(ss_2, effect_allele)))==1 & nchar(as.character(pull(ss_2, noneffect_allele)))==1),]

ss_2 <- ss_2[!((as.character(pull(ss_2, effect_allele))=="A" & as.character(pull(ss_2, noneffect_allele))=="T") | 
             (as.character(pull(ss_2, effect_allele))=="T" & as.character(pull(ss_2, noneffect_allele))=="A") | 
             (as.character(pull(ss_2, effect_allele))=="C" & as.character(pull(ss_2, noneffect_allele))=="G") |
             (as.character(pull(ss_2, effect_allele))=="G" & as.character(pull(ss_2, noneffect_allele))=="C")),]

duplicates <- ss_2[duplicated(pull(ss_2, rsID_col)), 1]
ss_2 <- ss_2[!(pull(ss_2, rsID_col)) %in% duplicates, ]

PRS_CS <- ss_2[, c(rsID_col, effect_allele, noneffect_allele, "beta_num", "P_num")]
colnames(PRS_CS) <- c("SNP", "A1", "A2", "BETA", "P")

write.table(ss_2, paste0(sum.stat, "_edited_by_pipeline_for_LDpred2"), row.names = F, quote = F, sep="\t")

write.table(PRS_CS, paste0(sum.stat, "_edited_by_pipeline_for_PRS_CS"), row.names = FALSE, quote = FALSE, sep="\t")
cat("Edited GWAS for PRS-CS and LDpred2....\nEditing GWAS completed!\n")

