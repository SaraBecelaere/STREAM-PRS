#!/usr/bin/bash

PRScs_path=$1
GWAS_file=$2
test_data=$3
test_prefix=$4
training_data=$5
training_prefix=$6
out_PRScs=$7
PRS_CS_ref_files=$8
GWAS_size=$9
phi_values=${10}

#Get number of cores
export MKL_NUM_THREADS=$SLURM_CPUS_ON_NODE
export NUMEXPR_NUM_THREADS=$SLURM_CPUS_ON_NODE
export OMP_NUM_THREADS=$SLURM_CPUS_ON_NODE

#Transform phi values
IFS=',' read -ra PHI_array <<< "$phi_values"

mkdir -p ${out_PRScs}/temp

#Loop through phi values for training and test

for phi_value in "${PHI_array[@]}"

do
	seq 1 22 | parallel -j 22 "python ${PRScs_path}/PRScs.py --ref_dir=${PRS_CS_ref_files} --bim_prefix=${training_data} --sst_file=${GWAS_file} --n_gwas=${GWAS_size} --chrom={} --phi=${phi_value}  --out_dir=${out_PRScs}/temp/${training_prefix}"
	
	check_missing_files() {
    missing_files=()
    for i in {1..22}; do
        file="${out_PRScs}/temp/${training_prefix}_pst_eff_a1_b0.5_phi${phi_value}_chr${i}.txt"
        if [ ! -f "$file" ]; then
            missing_files+=("$i")
        fi
    done
	}
	
	j=1
	while true; do
	check_missing_files
	
		if [ ${#missing_files[@]} -eq 0 ]; then
			echo "PRS-CS completed for phi" ${phi_value}
			break
		
		elif [ $j -eq 10]; then
			echo "Tried PRS-CS ten times and still no estimates for all chromosomes... exiting..."
			break
			
		else
			# Output the missing files
			echo "The following chromosomes are missing for phi ${phi_value} ${missing_files[@]}"

			# Run the Python script for each missing file
			for number in "${missing_files[@]}"; do
				echo "Running PRScs.py for chromosome ${number}..."
				python "${PRScs_path}/PRScs.py" --ref_dir="${PRS_CS_ref_files}" --bim_prefix="${training_data}" --sst_file="${GWAS_file}" --n_gwas="${GWAS_size}" --chrom="${number}" --phi="${phi_value}" --out_dir="${out_PRScs}/temp/${training_prefix}"
			done
			((j++))
		fi
	done

	for i in {1..22}
	do
		cat ${out_PRScs}/temp/${training_prefix}_pst_eff_a1_b0.5_phi${phi_value}_chr${i}.txt >> ${out_PRScs}/${training_prefix}_pst_eff_a1_b0.5_phi${phi_value}_all.txt
	done

	plink --bfile ${training_data} --score ${out_PRScs}/${training_prefix}_pst_eff_a1_b0.5_phi${phi_value}_all.txt 2 4 6 --out ${out_PRScs}/${training_prefix}_phi_${phi_value}
	plink --bfile ${test_data} --score ${out_PRScs}/${training_prefix}_pst_eff_a1_b0.5_phi${phi_value}_all.txt 2 4 6 --out ${out_PRScs}/${test_prefix}_phi_${phi_value}

done

rm -r ${out_PRScs}/temp
