
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


# at the model fit


# Simulate data
#--------------
DIR=/hpc/users/hoffmg01/work/lucida_analysis/simulations
# FIT=$DIR/fit_lucida.RDS
# DATA=$DIR/lucida_fit_data.RDS
FIT=$DIR/test_lucida_fit.RDS
DATA=$DIR/test_lucida_fit_data.RDS

echo "" > script_sim.sh
for libScaleFactor in $(echo 0.5 1 5 10 25)
do
for i in $(seq 1 1 25)
do
  ID=${libScaleFactor}_${i}
  echo "$DIR/create_dataset.R --fit $FIT --data $DATA --subject donor_id --seed $i --logFC 0.07 --pDE 0.05 --libScaleFactor ${libScaleFactor} --output /sc/arion/scratch/hoffmg01/sims/sim_${ID}.h5ad" >> script_sim.sh
done
done

cat script_sim.sh | parallel -P 40


# Recode for efficient access
#----------------------------
SRC=/hpc/users/hoffmg01/work/GenomicDataStream_analysis/recode_h5ad.py 

echo "" > script_recode.sh
for libScaleFactor in $(echo 0.5 1 5 10 25)
do
for i in $(seq 1 1 25)
do
  ID=${libScaleFactor}_${i}
  FILE=/sc/arion/scratch/hoffmg01/sims/sim_${ID}.h5ad
  OUT=/sc/arion/scratch/hoffmg01/sims/sim_${ID}_recode.h5ad
  echo "$SRC --input $FILE --sortBy cell_type,donor_id --format CSC --compression lzf --out $OUT " >> script_recode.sh
done
done

cat script_recode.sh | parallel -P 40

# Run DE analysis
#----------------

echo "" > script_de.sh
for libScaleFactor in $(echo 0.5 1 5 10 25)
do
for i in $(seq 1 1 25)
do
  ID=${libScaleFactor}_${i}
  FILE=/sc/arion/scratch/hoffmg01/sims/sim_${ID}_recode.h5ad
  OUT=/sc/arion/scratch/hoffmg01/sims/res_sim_${ID}.parquet
  echo "$DIR/run_analysis.R --h5ad $FILE --formula \"~ Dx + (1|donor_id)\" --cluster_id cell_type --output $OUT" >>  script_de.sh
done
done

cat script_de.sh | parallel -P 40

# Performance plots
###################

rmarkdown::render("plot_results.Rmd")


methods <- c(  
    "lucida",
    "lucida [1 step]",
    "lucida [pb]",
    "glmGamPoi")

