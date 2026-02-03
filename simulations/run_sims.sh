
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

CTs = c("gamma-delta T cell", "mucosal invariant T cell", "naive B cell", "naive thymus-derived CD4-positive, alpha-beta T cell")

n_donors_array = c(4, 10, 25, 50, 100, 250, 400, 500, 700, 981)

for(n_donors in n_donors_array[c(1:5)]){

  message(n_donors)
  sce2 = sce[,sce$cell_type %in% CTs]

  if( n_donors < 25 ){
    sce2 = sce2[,sce2$sex == "female"]
  }

  sce2 = sce2[,sce2$donor_id %in% levels(sce2$donor_id)[seq(n_donors)]]
  colData(sce2) = droplevels(colData(sce2))

  fit = lucida(sce2, ~ age + sex + (1|donor_id), cluster_id = "cell_type")

  file = paste0("test_lucida_fit_", n_donors, ".RDS")
  saveRDS(fit, file=file)

  file = paste0("test_lucida_fit_data_", n_donors, ".RDS")
  saveRDS(colData(sce2), file=file)
}



# Simulate data
#--------------
ml parallel
DIR=/hpc/users/hoffmg01/work/lucida_analysis/simulations

NREPS=5
echo "" > script_sim.sh
for N in $(echo "4 10 25 50 100 250 400 500 700 981")
do
for libScaleFactor in $(echo 0.5 1 5 10 25)
do
for i in $(seq 1 1 $NREPS)
do
  FIT=$DIR/test_lucida_fit_${N}.RDS
  DATA=$DIR/test_lucida_fit_data_${N}.RDS
  ID=${N}_${libScaleFactor}_${i}
  echo "$DIR/create_dataset.R --fit $FIT --data $DATA --subject donor_id --seed $i --logFC 0.07 --pDE 0.05 --libScaleFactor ${libScaleFactor} --output /sc/arion/scratch/hoffmg01/sims/sim_${ID}.h5ad" >> script_sim.sh
done
done
done

cat script_sim.sh | parallel -P 20


# Recode for efficient access
#----------------------------
SRC=/hpc/users/hoffmg01/work/GenomicDataStream_analysis/recode_h5ad.py 

echo "" > script_recode.sh
for N in $(echo "4 10 25 50 100 250 400 500 700 981")
do
for libScaleFactor in $(echo 0.5 1 5 10 25)
do
for i in $(seq 1 1 $NREPS)
do
  ID=${N}_${libScaleFactor}_${i}
  FILE=/sc/arion/scratch/hoffmg01/sims/sim_${ID}.h5ad
  OUT=/sc/arion/scratch/hoffmg01/sims/sim_${ID}_recode.h5ad
  echo "$SRC --input $FILE --sortBy cell_type,donor_id --format CSC --compression lzf --out $OUT " >> script_recode.sh
done
done
done

cat script_recode.sh | parallel -P 40

# Run DE analysis
#----------------

echo "" > script_de.sh
for N in $(echo "4 10 25 50 100 250 400 500 700 981")
do
for libScaleFactor in $(echo 0.5 1 5 10 25)
do
for i in $(seq 1 1 $NREPS)
do
  ID=${N}_${libScaleFactor}_${i}
  FILE=/sc/arion/scratch/hoffmg01/sims/sim_${ID}_recode.h5ad
  OUT=/sc/arion/scratch/hoffmg01/sims/res_sim_${ID}.parquet
  echo "$DIR/run_analysis.R --h5ad $FILE --formula \"~ Dx + (1|donor_id)\" --cluster_id cell_type --output $OUT" >>  script_de.sh
done
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

