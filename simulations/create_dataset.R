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

library(lucida)
library(SingleCellExperiment)
library(tidyverse)

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

library(anndataR)

write_h5ad(sce.sim, opt$output)