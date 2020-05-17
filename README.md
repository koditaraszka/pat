# PAT: pleiotropic association test
GWAS method for joint analysis of multiple traits

To use PAT follow these steps:
1. `git clone https://github.com/koditaraszka/pat`
2. `cd pat`
3. `python3 main.py [arguments]`

To see the list of arguments `python3 main.py -h`

## Required Arguments:
`-e/--envir`: text file where each line is whitespace separated name of two traits and their pairwise environmental correlation
`-g/--genetic`: text file where each line is whitespace separated name of two traits, genetic variance trait 1 and trait 2, and their genetic correlation
`-n/--pop_sizes`: text file where each line is the names of two traits, sample size for trait 1 and trait 2, and the overlapping individuals
`-f/--gwas_files`: text file where each line is the name of a trait and the path to its summary statistics

## Optional Arguments:
`-s/--sampling`: Scaling factor r for use in null simulations. Default is 8
`-m/--sims`: flag to simulate data rather than analyze --gwas_files
`-x/--num`: number of simulations to generate. Default: 1e6 (only set when `-m/--sims` is used)
`-y/--polygenic`: number of snps under polygenic model. If not passed number of snps in real data or simulations is used.
`-u/--null`: number of null simulations to run. Default: 1e6
`-z/--set-zeros`: each line is the name of a trait whose genetic effect to be set to zero in simulations (used only when `-m/--sims` is used)
`-i/--migwas`: flag indicating MI GWAS results should be computed as well
`-w/--sig_thresh`: Significance threshould for M-values. Default: 5e-8

# Directories:
The data and scripts should be readable but may be a little messy. If there are any questions please reach out.

## input:
this contains example input files and a readme.txt

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


