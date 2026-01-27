
# Lucida fit
#-----------

library(lucida)
library(SingleCellExperiment)

file = "/sc/arion/projects/CommonMind/hoffman/scRNAseq_data/Yazar_Science_2022/1k1k_sort_CSC_lzf.h5ad"
sce = readH5AD(file)
sce$libSize = sce$nCount_RNA

# filter to include only cell types with 1k observations
tab = table(sce$cell_type)
sce = sce[,sce$cell_type %in% names(tab[tab>1000])]
sce = sce[,sce$libSize > 1000]
sce$cell_type = droplevels(sce$cell_type)

CT = levels(sce$cell_type)[1:2]
sce2 = sce[,sce$cell_type %in% CTs]
sce2 = sce2[,sce2$donor_id %in% levels(sce2$donor_id)[1:300]]

fit = lucida(sce2, ~ age + sex + (1|donor_id), 
  cluster_id = "cell_type")

saveRDS(fit, file="test_lucida_fit.RDS")
saveRDS(colData(sce2), file="test_lucida_fit_data.RDS")



# fit2 = lucida(sce.sim, ~ Dx + (1|donor_id), cluster_id = "cell_type")



# Simulate data
#--------------
DIR=/hpc/users/hoffmg01/work/lucida_analysis/simulations
# FIT=$DIR/fit_lucida.RDS
# DATA=$DIR/lucida_fit_data.RDS
FIT=$DIR/test_lucida_fit.RDS
DATA=$DIR/test_lucida_fit_data.RDS

$DIR/create_dataset.R --fit $FIT --data $DATA --subject donor_id --seed 1 --logFC 0.07 --pDE 0.3 --libScaleFactor 1 --output /sc/arion/scratch/hoffmg01/sims/sim_1.h5ad


# Recode for efficient access
#----------------------------
SRC=/hpc/users/hoffmg01/work/GenomicDataStream_analysis/recode_h5ad.py 
FILE=/sc/arion/scratch/hoffmg01/sims/sim_1.h5ad

$SRC --input $FILE --sortBy cell_type,donor_id --format CSC --compression lzf --out /sc/arion/scratch/hoffmg01/sims/sim_1_recode.h5ad


# Run DE analysis
#----------------
run_analysis.R --h5ad /sc/arion/scratch/hoffmg01/sims/sim_1_recode.h5ad --formula "~ Dx + (1|donor_id)" --cluster_id cell_type --output /sc/arion/scratch/hoffmg01/sims/res_sim_1.parquet


1) DE for somatic rate
2) check anndataR written h5ad, can be read in python?
  uns versus obs