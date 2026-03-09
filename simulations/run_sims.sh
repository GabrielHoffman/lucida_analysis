
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

CTs = c("mucosal invariant T cell", "gamma-delta T cell", "naive B cell", "naive thymus-derived CD4-positive, alpha-beta T cell")

n_donors_array = c(10, 25, 50, 100, 250, 400, 500, 700, 981)

for(n_donors in n_donors_array){

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

# testing
NREPS=10
NSAMPLES="10 25 50"  
LSF="1" # libScaleFactors
OUTFOLDER=/sc/arion/scratch/hoffmg01/sims/1k1k_v1/
LOGFC=0.2

# Production
# NREPS=50
# NSAMPLES="10 25 50 100 250 400 500"  # 700 981
# LSF="0.5 1 5 10 25" # libScaleFactors
# OUTFOLDER=/sc/arion/scratch/hoffmg01/sims/1k1k_v1/
# LOGFC=0.07 

mkdir -p $OUTFOLDER

echo "" > script_sim.sh
for N in $(echo $NSAMPLES)
do
for libScaleFactor in $(echo $LSF)
do
for i in $(seq 1 1 $NREPS)
do
  FIT=$DIR/test_lucida_fit_${N}.RDS
  DATA=$DIR/test_lucida_fit_data_${N}.RDS
  ID=${N}_${libScaleFactor}_${i}
  echo "$DIR/create_dataset.R --fit $FIT --data $DATA --subject donor_id --seed $i --logFC $LOGFC --pDE 0.05 --libScaleFactor ${libScaleFactor} --output $OUTFOLDER/sim_${ID}.h5ad" >> script_sim.sh
done
done
done

# run sims
cat script_sim.sh | parallel -P 60

# check that files were written
cat script_sim.sh | awk '{print $NF}' | xargs ls > /dev/null


# Recode for efficient access
#----------------------------
SRC=/hpc/users/hoffmg01/work/GenomicDataStream_analysis/recode_h5ad.py 

echo "" > script_recode.sh
for N in $(echo $NSAMPLES)
do
for libScaleFactor in $(echo $LSF)
do
for i in $(seq 1 1 $NREPS)
do
  ID=${N}_${libScaleFactor}_${i}
  FILE=$OUTFOLDER/sim_${ID}.h5ad
  OUT=$OUTFOLDER/sim_${ID}_recode.h5ad
  echo "$SRC --input $FILE --sortBy cell_type,donor_id --format CSC --compression lzf --out $OUT " >> script_recode.sh
done
done
done

cat script_recode.sh | parallel -P 25

# check that files were written
cat script_recode.sh | awk '{print $NF}' | xargs ls > /dev/null

# remove tmp h5ad's
comm -3 <(ls $OUTFOLDER/*recode.h5ad | sort) <(ls $OUTFOLDER/*.h5ad | sort) | xargs -n1 rm -f




# Run DE analysis
#----------------

# rm -f /sc/arion/scratch/hoffmg01/sims/1k1k_v1/*.parquet

METHODS=/hpc/users/hoffmg01/work/lucida_analysis/simulations/methods.in

echo "" > script_de.sh
for N in $(echo $NSAMPLES)
do
for libScaleFactor in $(echo $LSF)
do
for i in $(seq 1 1 $NREPS)
do
  ID=${N}_${libScaleFactor}_${i}
  FILE=$OUTFOLDER/sim_${ID}_recode.h5ad
  OUT=$OUTFOLDER/res_sim_${ID}.parquet
  echo "$DIR/run_analysis.R --h5ad $FILE --formula \"~ Dx + (1|donor_id)\" --cluster_id cell_type --methods $METHODS --output $OUT" >>  script_de.sh
done
done
done

cat script_de.sh | parallel -P 60

# check that files were written
cat script_de.sh | awk '{print $NF}' | xargs ls > /dev/null 2> err.log

cat err.log | tr "'" " " | awk '{print $4}' | parallel basename {} .parquet > jobs.prefix

grep -f jobs.prefix script_de.sh | parallel -P60


cat script_de.sh | grep -v "sim_25_\|sim_50_\|sim_100_\|sim_250_" | parallel -P60


cat script_de.sh | sed 's/res_sim_/test2\/res_sim_/g' | parallel



# Performance plots
###################

rmarkdown::render("plot_results.Rmd")










# TRAJECTORY

# Simulate data
#--------------
ml parallel
DIR=/hpc/users/hoffmg01/work/lucida_analysis/simulations

# testing
NREPS=10
NSAMPLES="10 25 50"  
LSF="1" # libScaleFactors
OUTFOLDER=/sc/arion/scratch/hoffmg01/sims/1k1k_v1/trajectory
LOGFC=0.2



# Production
# NREPS=50
# NSAMPLES="10 25 50 100 250 400 500"  # 700 981
# LSF="0.5 1 5 10 25" # libScaleFactors
# OUTFOLDER=/sc/arion/scratch/hoffmg01/sims/1k1k_v1/
# LOGFC=0.07 

mkdir -p $OUTFOLDER

echo "" > script_sim.sh
for N in $(echo $NSAMPLES)
do
for libScaleFactor in $(echo $LSF)
do
for i in $(seq 1 1 $NREPS)
do
  FIT=$DIR/test_lucida_fit_${N}.RDS
  DATA=$DIR/test_lucida_fit_data_${N}.RDS
  ID=${N}_${libScaleFactor}_${i}
  echo "$DIR/create_dataset.R --fit $FIT --data $DATA --subject donor_id --seed $i --logFC $LOGFC --pDE 0.05 --libScaleFactor ${libScaleFactor} --output $OUTFOLDER/sim_${ID}.h5ad" >> script_sim.sh
done
done
done

# run sims
cat script_sim.sh | parallel -P 60

# check that files were written
cat script_sim.sh | awk '{print $NF}' | xargs ls > /dev/null


# Recode for efficient access
#----------------------------
SRC=/hpc/users/hoffmg01/work/GenomicDataStream_analysis/recode_h5ad.py 

echo "" > script_recode.sh
for N in $(echo $NSAMPLES)
do
for libScaleFactor in $(echo $LSF)
do
for i in $(seq 1 1 $NREPS)
do
  ID=${N}_${libScaleFactor}_${i}
  FILE=$OUTFOLDER/sim_${ID}.h5ad
  OUT=$OUTFOLDER/sim_${ID}_recode.h5ad
  echo "$SRC --input $FILE --sortBy cell_type,donor_id --format CSC --compression lzf --out $OUT " >> script_recode.sh
done
done
done

cat script_recode.sh | parallel -P 25

# check that files were written
cat script_recode.sh | awk '{print $NF}' | xargs ls > /dev/null

# remove tmp h5ad's
comm -3 <(ls $OUTFOLDER/*recode.h5ad | sort) <(ls $OUTFOLDER/*.h5ad | sort) | xargs -n1 rm -f




# Run DE analysis
#----------------

# rm -f /sc/arion/scratch/hoffmg01/sims/1k1k_v1/*.parquet

echo "" > script_de.sh
for N in $(echo $NSAMPLES)
do
for libScaleFactor in $(echo $LSF)
do
for i in $(seq 1 1 $NREPS)
do
  ID=${N}_${libScaleFactor}_${i}
  FILE=$OUTFOLDER/sim_${ID}_recode.h5ad
  OUT=$OUTFOLDER/res_sim_${ID}.parquet
  echo "$DIR/run_analysis.R --h5ad $FILE --formula \"~ Dx + (1|donor_id)\" --cluster_id cell_type --output $OUT" >>  script_de.sh
done
done
done

cat script_de.sh | parallel -P 60

# check that files were written
cat script_de.sh | awk '{print $NF}' | xargs ls > /dev/null 2> err.log

cat err.log | tr "'" " " | awk '{print $4}' | parallel basename {} .parquet > jobs.prefix

grep -f jobs.prefix script_de.sh | parallel -P60


cat script_de.sh | grep -v "sim_25_\|sim_50_\|sim_100_\|sim_250_" | parallel -P60


cat script_de.sh | sed 's/res_sim_/test2\/res_sim_/g' | parallel



# Performance plots
###################

rmarkdown::render("plot_results.Rmd")











