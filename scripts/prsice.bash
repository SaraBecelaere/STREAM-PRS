#!/usr/bin/bash

PRSice_path=$1
GWAS_file=$2
training_data=$3
training_prefix=$4
snp_col=$5
chromosome_col=$6
pos_col=$7
effect_allele=$8
noneffect_allele=$9
clumpkb=${10}
clumpp=${11}
clumpr2=${12}
pval_thresholds=${13}
output_PRSice=${14}
test_data=${15}
test_prefix=${16}

echo "Start PRSice-2..."

#PRSice for training
Rscript ${PRSice_path}/PRSice.R --prsice ${PRSice_path}/PRSice_linux --base ${GWAS_file} --target ${training_data} --snp ${snp_col} --chr ${chromosome_col} --bp ${pos_col} --A1 ${effect_allele} --A2 ${noneffect_allele} --stat beta_num --beta --pvalue P_num --clump-kb ${clumpkb} --clump-p ${clumpp} --clump-r2 ${clumpr2} --bar-levels ${pval_thresholds} --fastscore --no-regress --print-snp --out ${output_PRSice}/${training_prefix}

echo "PRSice-2 for training data completed...."

#PRSice for test
Rscript ${PRSice_path}/PRSice.R --prsice ${PRSice_path}/PRSice_linux --base ${GWAS_file} --target ${test_data} --snp ${snp_col} --chr ${chromosome_col} --bp ${pos_col} --A1 ${effect_allele} --A2 ${noneffect_allele} --stat beta_num --beta --pvalue P_num --extract ${output_PRSice}/${training_prefix}.snp --no-clump --bar-levels ${pval_thresholds} --fastscore --no-regress --print-snp --out ${output_PRSice}/${test_prefix}

echo "PRSice-2 for test data completed...."
echo "PRSice-2 completed!"
