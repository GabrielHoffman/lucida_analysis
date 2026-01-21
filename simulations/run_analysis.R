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
source("de_methods.R")

# read data
sce.sim <- readH5AD( opt$h5ad )

# select methods to run
methods <- c(  
    "lucida",
    "lucida [1 step]",
    "lucida [Bayesian]",
    "lucida [pb]",
    "nebula",
    "nebula (HL)" ,
    "dreamlet" ,
    "DESeq2" ,
    "edgeR",
    "glmGamPoi",
    "MAST")

# Run analysis
df <- run_analysis(sce.sim, 
      formula = as.formula(opt$formula), 
      cluster_id = opt$cluster_id, 
      methods = methods[-c(5,6)])

# write results to parquet
library(arrow)
write_parquet( df, sink = opt$out)
