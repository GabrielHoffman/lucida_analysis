
# Lucida fit
#-----------

cd /hpc/users/hoffmg01/work/lucida_analysis/simulations

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

CTs = c("gamma-delta T cell", "mucosal invariant T cell")
sce2 = sce[,sce$cell_type %in% CTs]
sce2 = sce2[,sce2$donor_id %in% levels(sce2$donor_id)[1:300]]

fit = lucida(sce2, ~ age + sex + (1|donor_id), 
  cluster_id = "cell_type")

saveRDS(fit, file="test_lucida_fit.RDS")
saveRDS(colData(sce2), file="test_lucida_fit_data.RDS")





# Simulate data
#--------------
DIR=/hpc/users/hoffmg01/work/lucida_analysis/simulations
# FIT=$DIR/fit_lucida.RDS
# DATA=$DIR/lucida_fit_data.RDS
FIT=$DIR/test_lucida_fit.RDS
DATA=$DIR/test_lucida_fit_data.RDS
 
for i in $(seq 1 1 10)
do
  $DIR/create_dataset.R --fit $FIT --data $DATA --subject donor_id --seed $i --logFC 0.07 --pDE 0.1 --libScaleFactor 1 --output /sc/arion/scratch/hoffmg01/sims/sim_${i}.h5ad &
done

# Recode for efficient access
#----------------------------
SRC=/hpc/users/hoffmg01/work/GenomicDataStream_analysis/recode_h5ad.py 

for i in $(seq 1 1 10)
do
  FILE=/sc/arion/scratch/hoffmg01/sims/sim_${i}.h5ad
  OUT=/sc/arion/scratch/hoffmg01/sims/sim_${i}_recode.h5ad
  $SRC --input $FILE --sortBy cell_type,donor_id --format CSC --compression lzf --out $OUT
done

# Run DE analysis
#----------------
for i in $(seq 1 1 10)
do
  FILE=/sc/arion/scratch/hoffmg01/sims/sim_${i}_recode.h5ad
  OUT=/sc/arion/scratch/hoffmg01/sims/res_sim_${i}.parquet
  $DIR/run_analysis.R --h5ad $FILE --formula "~ Dx + (1|donor_id)" --cluster_id cell_type --output $OUT
done

# Performance plots
###################

library(arrow)
library(tidyverse)

file = "/sc/arion/scratch/hoffmg01/sims/res_sim_1.parquet"
df = read_parquet(file)





