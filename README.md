# STREAM-PRS
STREAM-PRS (**S**treamlined **T**oolkit for **R**eliable **E**valuation and **A**nalysis of **M**ultiple **P**olygenic **R**isk **S**cores) is a streamlined pipeline for calculating polygenic risk scores (PRS). It performs:
1. Basic editing of GWAS summary statistics
2. Calculation of PRS with the following tools: PRSice-2, PRS-CS, LDpred2, lassosum and lassosum2.
3. PC-correction and standardization of the obtained PRS for all tools.
4. Selection of the best PRS per tool and overall

## Getting started
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

## How to use

To use the PRS pipeline, you need to edit the STREAM-PRS.bash file. Below you will find a list of all the flags you need to fill in, and which are optional.

### Essential parameters



## Citation

## Contact

