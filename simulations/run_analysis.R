#! /usr/bin/env Rscript

library(getopt)

spec <- matrix(c(
  'h5ad',       'i', 1, "character",
  'formula',    'f', 1, "character",
  'cluster_id', 'c', 1, "character",
  'output',     'o', 1, "character"
), byrow=TRUE, ncol=4)

# Parse command line arguments
opt <- getopt(spec)

# load analysis code
message("Loading packages...")
source("/hpc/users/hoffmg01/work/lucida_analysis/simulations/de_methods.R")

# read data
message("Read H5AD...")
sce.sim <- readH5AD( opt$h5ad )

# select methods to run
methods <- c(  
    "lucida",
    "lucida [1 step]",
    "lucida [Bayesian]",
    "lucida [pb]",
    # "nebula",
    # "nebula (HL)" ,
    "dreamlet" ,
    "DESeq2" ,
    "edgeR",
    "glmGamPoi",
    "MAST")

# Run analysis
message("Run analysis...")
df <- run_analysis(sce.sim, 
      formula = as.formula(opt$formula), 
      cluster_id = opt$cluster_id, 
      methods = methods)

# Get sim params
file <- gsub("_recode.h5ad$", ".info.RDS", opt$h5ad)
params <- readRDS(file)

# identifier for this simulation
simID <- paste(params, collapse="_")

# write results to parquet
library(arrow)

df %>%
  mutate(simID = simID,
    cluster_id = factor(cluster_id),
    ID = factor(ID)) %>%
  write_parquet( sink = opt$out)
