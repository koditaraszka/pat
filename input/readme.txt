Input files for PAT

Each line of the input files needs the trait names.
Trait names are dictionary keys and can be set by the user.
The key used will be part of the output and must be consistent across input files.

Required:

env.txt (passed to -e/--envir) contains the environmental correlation between traits
The value here is derived by weighting the cross-trait LD-Score intercept by sample overlap and subtracting out the genetic covariance
Do not reweight the environmental correlation by sample overlap, this is done internally
Each line should be space delimited with: trait1 trait2 envir_corr
The order of the trait names does not matter as long as all pairs of traits are present

gen.txt (passed to -g/--genetic) contains the genetic variance and covariance between traits
These values come directly from LD-Score and cross-trait LD-Score regression
Each line should be space delimited with: trait1 trait2 gen_var1 gen_var2 gen_corr
Again the order does not matter as long as all pairs of traits are present

file.txt (passed to -f/--gwas_files) matches the trait name (key) to the input file
Each line should be space delimited with: trait input_file.txt
The input files are assumed to have these columns with these names: RSID, CHR, BP, A1, A2, Z, P
Additional rows are allowed and will be given appropriate suffixes to indicate the source file.
For example if you have the column N, the output will have the column N_trait

pop.txt (passed to -n/--pop_sizes) contains the population sample sizes
Each line should be space delimited with: trait1 trait2 sample_size1 sample_size2 overlap
We use the smaller sample size to estimate overlap, but if the value is a user input

Optional:

z_*.txt (passed to -z/--set-zeros) which lists the traits to set genetic effect to zero
Each line is the name of a trait (Again this must be the same key across the files)
This argument is only used in simulations (when -i/--sims is called)
