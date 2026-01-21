


# Simulate data
DIR=/hpc/users/hoffmg01/work/lucida_analysis/simulations
FIT=$DIR/fit_lucida.RDS
DATA=~/work/lucida_analysis/simulations/lucida_fit_data.RDS


# $DIR/create_dataset.R

R --args --fit $FIT --data $DATA --subject donor_id --seed 1 --logFC 0.07 --pDE 0.3 --libScaleFactor 1 --output /sc/arion/scratch/hoffmg01/sims/sim_1.h5ad


# Recode for efficient access
#----------------------------
SRC=/hpc/users/hoffmg01/work/GenomicDataStream_analysis/recode_h5ad.py 
FILE=/sc/arion/scratch/hoffmg01/sims/sim_1.h5ad

$SRC --input $FILE --sortBy cluster_id,donor_id --format CSC --compression lzf --out /sc/arion/scratch/hoffmg01/sims/sim_1_recode.h5ad


# Run DE analysis
#----------------
run_analysis.R --h5ad /sc/arion/scratch/hoffmg01/sims/sim_1_recode.h5ad --formula "~ Dx + (1|donor_id)" --cluster_id "cluster_id" --output /sc/arion/scratch/hoffmg01/sims/res_sim_1.parquet