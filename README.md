# STREAM-PRS
STREAM-PRS (**S**treamlined **T**oolkit for **R**eliable **E**valuation and **A**nalysis of **M**ultiple **P**olygenic **R**isk **S**cores) is a streamlined pipeline for calculating polygenic risk scores (PRS). It performs:
1. Basic editing of GWAS summary statistics
2. Calculation of PRS with the following tools: PRSice-2, PRS-CS, LDpred2, lassosum and lassosum2.
3. PC-correction and standardization of the obtained PRS for all tools.
4. Selection of the best PRS per tool and overall

## Getting started

### Clone repository and install conda environment

First, clone the repository as follows:

```
git clone https://github.com/SaraBecelaere/STREAM-PRS.git
```

Next, we will install the conda environment. See [here](https://docs.anaconda.com/miniconda/install/) to find instructions on how to install miniconda if you haven't already. Go to the `conda_env` folder and install the conda environment.

```
cd conda_env
conda env create -f STREAM-PRS.yml
``` 

To activate the conda environment use:

```
conda activate STREAM-PRS
```

### Install tools and download necessary files

#### PRSice-2

- Tool: Download an executable for your system from this [website](https://choishingwan.github.io/PRSice/) and unzip the folder. E.g. for linux:

    ```
    unzip PRSice_linux.zip -d PRSice
    ```

- Reference files: not necessary to download for PRSice-2
  
#### PRS-CS

- Tool: Download the GitHub page as follows:

    ```
    git clone https://github.com/getian107/PRScs.git
    ```

- Reference files: From their GitHub page, download the reference data that is similar to your dataset. E.g. for European data: click [here](https://www.dropbox.com/scl/fi/zcq9ocvzp2bcfg05ug7lg/ldblk_ukbb_eur.tar.gz?rlkey=e3eaetxh4b1faennmom4e3mmi&e=1&dl=0) and download the data, then use the following code:

  ```
  tar -zxvf ldblk_ukbb_eur.tar.gz
  ```

  For other reference populations, follow the instructions on the [PRS-CS GitHub page](https://github.com/getian107/PRScs).

#### LDpred2 and lassosum2

- Tool: not necessary to download, they are integrated in the STREAM-PRS conda environment
  
- Reference files: Use the following [link]([https://figshare.com/articles/dataset/LD_reference_for_HapMap3_/21305061](https://drive.usercontent.google.com/download?id=17dyKGA2PZjMsivlYb_AjDmuZjM1RIGvs&export=download&authuser=0)) to download the LD reference matrices. Download `ldref_hm3_plus.zip`. Then unzip the file:

  ```
  unzip ldref_hm3_plus.zip -d ldref_hm3_plus
  ```

#### lassosum

- Tool: not necessary to download, the tool is integrated in the STREAM-PRS conda environment

- Reference files: not necessary to download for lassosum

## How to use

To use the PRS pipeline, you need to edit the STREAM-PRS.bash file. Below you will find a list of all the flags you need to fill in, and which are optional.

### Parameters that you have to adapt

#### Tools

- PRSice_path: fill in the full path to the PRSice folder, where the PRSice.R and (for linux) PRSice_linux are stored
- PRScs_path: fill in the full path to the PRScs folder, where PRScs.py is stored

#### Model

- binary_trait: fill in either TRUE (the trait you calculate PRS for is binary) or FALSE (the trait you calculate PRS for is quantitative)
- cross-validation: fill in either TRUE (you want to use cross-validation) or FALSE (you do not want to use cross-validation)

#### Output path

- out: fill in the full output path here (. if you want to use the current folder)

#### GWAS summary statistics

- GWAS_file: fill in the full path to the GWAS summary statistics file you want to use, where the variant IDs are stored as chr_pos_ID
- GWAS_file_rsID: fill in the full path to the GWAS summary statistics file you want to use, where the variant IDs are stored as rsID

*Note: if the GWAS summary statistics file you use for chr_pos_ID and rsID are the same, use the same file for GWAS_file and GWAS_file_rsID*

#### GWAS column names

- snp_col: fill in column name that indicates the variant ID in chr_pos format in the GWAS summary statistics file
- rsID_col: fill in column name that indicates the variant ID in rsID format in the GWAS summary statistics file
- chromosome_col: fill in column name that indicates the chromosome in the GWAS summary statistics file
- pos_col: fill in column name that indicates the basepair position in the GWAS summary statistics file
- effect_allele: fill in column name that indicates the effect allele in the GWAS summary statistics file
- noneffect_allele: fill in column name that indicates the non effect allele in the GWAS summary statistics file
- beta_col: fill in column name that indicates the SNP effect size in beta format in the GWAS summary statistics file
- beta_se_col: fill in column name that indicates the standard error of the beta in the GWAS summary statistics file
- p_value_col: fill in column name that indicates the P-value of the SNP in the GWAS summary statistics file
- N_col: fill in column name that indicates the sample size of the SNP in the GWAS summary statistics file (*Note: if this column is not provided in the GWAS summary statistics, add this column yourself. Just use the GWAS size for all SNPs.*)
- allele_freq_col: fill in column name that indicates the allele frequency of the effect allele in the GWAS summary statistics file (*Note: this is only necessary of you use the alternative LDpred2 script - see below - if you do not use this, just fill in "" instead*)

#### Size of the GWAS

- GWAS_size: fill in the GWAS size number (cases and controls combined) here e.g. 100000

#### Training data

*Note: training data should be in plink .bed, .bim, .fam format*

- training_file: fill in the full path to the training plink files in chr_pos variant ID format
- training_file_prefix: fill in plink file prefix (only prefix!!) of the training data in chr_pos variant ID format
- training_file_rsID: fill in the full path to the training plink files in rsID variant ID format
- training_file_rsID_prefix: fill in plink file prefix (only prefix!!) of the training data in rsID variant ID format

#### Test data

*Note: test data should be in plink .bed, .bim, .fam format*

- test_file: fill in the full path to the test plink files in chr_pos variant ID format
- test_file_prefix: fill in plink file prefix (only prefix!!) of the test data in chr_pos variant ID format
- test_file_rsID: fill in the full path to the test plink files in rsID variant ID format
- test_file_rsID_prefix: fill in plink file prefix (only prefix!!) of the test data in rsID variant ID format

#### Phenotype file

- pheno: fill in the full path to the phenotype file

*Note: the phenotype file should be in the following format: FID IID PHENO and shouldn't include any other columns*

#### PC files

- PC_training: fill in the full path to the file containing the first x principal components of the training data
- PC_test: fill in the full path to the file containing the first x principal components of the test data

*Note: this files should be in the following format: FID IID PC1 PC2 ... PCX. The number of PCs that you use is optional, but the number should be the same for the training and test data. The file should only include FID, IID and PCs, no additional columns*

#### Cores

- cores: fill in the number of cores you use to run the pipeline

#### PRSice specific parameters

*Note: all PRSice parameters are optional to adapt.*

- pval_thresholds: fill in a comma separated list of P-value thresholds to test (default: 5e-08,1e-05,0.0001,0.001,0.05,0.1,0.5,1).
- clumpkb: fill in the size of the window that should be used for clumping, in kb (default: 250).
- clumpp: fill in the p-values that should be considered while clumping (default: 1.000000).
- clumpr2: fill in the correlation coefficient between the SNPs that is tolerated for clumping (default: 0.100000).

#### Lassosum specific parameters

*Note: Make sure that you specify ld_blocks_lasso and that it matches your population of interest and the genome build of your data. lasso_thresholding_values and lambda_values can be adapted, but this is optional.*

- ld_blocks_lasso: fill in which LD blocks you want to use (default: EUR.hg38).
- lasso_thresholding_values: fill in a comma separated list of shrinkage values that should be tested (default: 0.1,0.2,0.4,0.5,0.7,0.9,1)
- lambda_values: fill in which lambda values should be tested (default: "exp(seq(log(0.001), log(0.1), length.out = 20))").

#### PRS-CS specific parameters

*Note: Make sure that you specify reference_files_PRScs, adapting phi_values is optional.*

- reference_files_PRScs: fill in the full path to the directory that contains the reference files for PRS-CS
- phi_values: fill in a comma separated list of phi values that should be tested (default: 1e+00,1e-02,1e-04,1e-06).

#### LDpred2 and lassosum2 specific parameters

*Note: Make sure that you specifiy ldref_hm3_plus, adapting the other parameters is optional.*

- ldref_hm3_plus: fill in the full path to the directory that contains the LD reference matrices for LDpred2 and lassosum2
- values_h2_grid: fill in a list of values that will be used to change the estimated heritability value in the grid model (default: "c(0.3, 0.7, 1, 1.4)").
- values_p_grid: fill in which values for the proportion of causal variants (p) should be tested in the grid model (default: "seq_log(1e-5, 1, length.out = 21)").
- initial_p_auto_grid: fill in which initial values for p should be tested in the auto grid model (default: "seq_log(1e-4, 0.2, length.out = 30)").
- delta: fill in which shrinkage values should be tested (default: "c(0.001, 0.01, 0.1, 1)").
- nlambda: fill in how many lambda values should be tested (default: 30).
- lambda_min_ratio: fill in what the minimum ratio between lambda values should be (default: 0.01).

### Run STREAM-PRS

#### How to run

You can run STREAM-PRS with the following code. Make sure that you adapted the STREAM-PRS.bash file according to your data and folder structures.

```
bash STREAM-PRS.bash
```

#### Structure of STREAM-PRS.bash script

The file starts with all the parameters that should be filled in (see above). Then, it runs the following scripts sequentially:

- edit_GWAS.r: performs basic editing of the GWAS summary statistics files and creates the correctly formatted files that will serve as input for STREAM-PRS
- prsice.bash: runs PRSice-2
- lassosum.r: runs lassosum
- PRS_CS.bash: runs PRS-CS
- LDpred2_and_lassosum2.r: runs LDpred2 and lassosum2
- get_best_PRS.r **OR** get_best_PRS_linear.r (depending on if the phenotype is indicated to be binary or not): PC-corrects and standardizes all scores. Then selects the best score per tool and across all tools.

*Note: if STREAM-PRS crashed midway one of the scripts, you can put a '#' in front of the scripts that already ran successfully. To run get_best_PRS.r or get_best_PRS_linear.r scores should be available for **ALL FIVE TOOLS**.

## Citation

If you use STREAM-PRS, please cite our paper. Becelaere S. et al. STREAM-PRS: A multi-tool pipeline for streamlining Polygenic Risk Score computation (2025). Currently, the paper is not yet published.

## Contact

For questions or problems when running STREAM-PRS, please leave a message on GitHub or send an email to sara.becelaere@kuleuven.be.
