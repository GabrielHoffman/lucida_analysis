#! /usr/bin/env Rscript

library(getopt)

spec <- matrix(c(
  'fit',        'i', 1, "character",
  'cluster_id', 'c', 1, "character",
  'logFC',      'l', 1, "numeric",
  'pDE',        'p', 1, "numeric",
  'libScaleFactor', 's', 1, "numeric",
  'output',     'o', 1, "character"
), byrow=TRUE, ncol=4)

library(lucida)
library(SingleCellExperiment)
library(tidyverse)

# load lucida model
fit <- readRDS( opt$fit )

# Simulate based on fit
sce.sim <- simulateCountData(
  fit = fit, 
  formula = ~ (1|Individual), 
  data = colData(sce), 
  target = "Dx", 
  logFC = opt$logFC, 
  pDE = opt$pDE, 
  libScaleFactor = opt$libScaleFactor)

library(anndataR)

write_h5ad(sce.sim, opt$output, compression = "lzf")