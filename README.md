# PAT: pleiotropic association test
GWAS method for joint analysis of multiple traits

To use PAT follow these steps:
1. `git clone https://github.com/koditaraszka/pat`
2. `cd pat`
3. `python3 main.py [arguments]`

To see the list of arguments `python3 main.py -h`

## Required Arguments:

## Optional Arguments:


# Directories:
The data and scripts should be readable but may be a little messy. If there are any questions please reach out.

## results:
Here are the results comparing methods shown in the paper

### simulations
This directory contains results from two sets of simulations
1. 700ksims_70kcausal:
PAT, HIPO, and MI GWAS are compared on 700k sets of simulated z-scores with 70,000 being causal. There are 7 different configurations of genetic effect explored with different causal effect sizes.

2. migwasVpat
PAT and MI GWAS are compared under different pleitropic models. There are simulations with and without environmental correlation.

Many of the scripts used will be found here or in the directory scripts (input.py in scripts) was used to simulate the data.

### realdata
This directory contains the real data results for 2018 UKBB summary statistics
1. mtag
results and scripts for MTAG

2. hipo
results and scripts for HIPO

3. pat
results and scripts for pat. Includes scripts and data for replication

The input files are found in input/ukbb2018.zip and some data cleaning scripts can be found in the scripts directory

### scripts
This directory contains various scripts used to analyze and clean the data. Some scripts can be found in other directories when only used there.


