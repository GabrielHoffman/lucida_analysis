#! /usr/bin/env Rscript

library(getopt)

spec <- matrix(c(
  'fit',        'i', 1, "character",
  'data',       'd', 1, "character",
  'subject',    'u', 1, "character",
  'logFC',      'l', 1, "numeric",
  'pDE',        'p', 1, "numeric",
  'libScaleFactor', 's', 1, "numeric",
  'seed',       'e', 1, "integer",
  'output',     'o', 1, "character"
), byrow=TRUE, ncol=4)

# Parse command line arguments
opt <- getopt(spec)

suppressPackageStartupMessages({
library(lucida)
library(SingleCellExperiment)
library(tidyverse)
library(anndataR)
library(RhpcBLASctl)
})

omp_set_num_threads(1)

# load lucida model
fit <- readRDS( opt$fit )
data <- readRDS( opt$data )

# set seed
set.seed(opt$seed)

# Set disease status for simulations
data$Dx = "Control"
Dx_lvls = sample(levels(data[[opt$subject]]), nlevels(data[[opt$subject]])/2)
data$Dx[data[[opt$subject]] %in% Dx_lvls] = "Disease"

formula = as.formula(paste("~(1|", opt$subject, ")"))

# Simulate based on fit
sce.sim <- simulateCountData(
  fit = fit, 
  formula = formula, 
  data = data, 
  target = "Dx", 
  logFC = opt$logFC, 
  pDE = opt$pDE, 
  libScaleFactor = opt$libScaleFactor)

# save simulation parameters
metadata(sce.sim)$params$seed <- opt$seed
metadata(sce.sim)$params$nDonor <- nlevels(sce.sim$donor_id)

# write params to separate file
params <- metadata(sce.sim)$params
file = gsub("h5ad$", "info.RDS", opt$output)
saveRDS(params, file = file)

# only data.frame is readable to python
metadata(sce.sim) <- metadata(sce.sim)$info

# write to file
write_h5ad(sce.sim, opt$output, mode="w", compression="none")