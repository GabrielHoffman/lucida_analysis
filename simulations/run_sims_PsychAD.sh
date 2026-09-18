
# Lucida fit
#-----------

cd /hpc/users/hoffmg01/work/lucida_analysis/simulations

library(lucida)
library(SingleCellExperiment)
library(GenomicDataStream)
library(parallel)

file = "/sc/arion/projects/CommonMind/hoffman/scRNAseq_data/psychAD/MSSM_2024-02-01_16_17_sort_CSC_lzf.h5ad"
sce = readH5AD(file)

# keep autosomal genes
keep = rowData(sce)$gene_chrom %in% 1:22
sce = sce[keep,]

# filter to include only cell types with 1k observations
tab = table(sce$subclass)
sce = sce[,sce$subclass %in% names(tab[tab>1000])]
sce = sce[,sce$libSize > 1000]
sce$subclass = droplevels(sce$subclass)

CTs = c("EN_L6_CT", "IN_SST", "Astro", "Oligo")

n_donors_array = c(4, 6, 8, 10, 12, 16, 20, 25, 50, 100, 250, 400, 500, 700, 1000)

res = mclapply(n_donors_array, function(n_donors){

  message(n_donors)
  sce2 = sce[,sce$subclass %in% CTs]

  sce2 = sce2[,sce2$SubID %in% levels(sce2$SubID)[seq(n_donors)]]
  colData(sce2) = droplevels(colData(sce2))

  form = ~ Age + Sex + PMI + (1|SubID)
  fit = lucida(sce2, form, cluster_id = "subclass", nReaders=4, nthreads=8)

  file = paste0("fits/PsychAD/test_lucida_fit_", n_donors, ".RDS")
  saveRDS(fit, file=file)

  file = paste0("fits/PsychAD/test_lucida_fit_data_", n_donors, ".RDS")
  saveRDS(colData(sce2), file=file)
}, mc.cores=12)


# Simulate data
#--------------
ml parallel
DIR=/hpc/users/hoffmg01/work/lucida_analysis/simulations/

# testing
NREPS=10
# NSAMPLES="25 50 100 250 400 500"  
# NSAMPLES="4 6 8 10 12 16 20" 
NSAMPLES="8 12 16 20 25 50" 
NMAX=1000
LSF="1" # libScaleFactors
# LOGFC=0.07 # large N
LOGFC=.3 # small N
# LOGFC=0
COVARIATES="'Age + Sex + PMI'"
OUTFOLDER=/sc/arion/scratch/hoffmg01/sims/PsychAD/constant/$(echo $NSAMPLES | tr ' ' '_')_${LOGFC}

# rm -f $OUTFOLDER/* $OUTFOLDER/logs/* $OUTFOLDER/jobs/*
mkdir -p $OUTFOLDER $OUTFOLDER/logs
cd $OUTFOLDER

echo "" > $OUTFOLDER/script_sim.sh
for N in $(echo $NSAMPLES)
do
for libScaleFactor in $(echo $LSF)
do
for i in $(seq 1 1 $NREPS)
do
  FIT=$DIR/fits/PsychAD/test_lucida_fit_${NMAX}.RDS
  DATA=$DIR/fits//PsychAD/test_lucida_fit_data_${N}.RDS
  ID=${N}_${libScaleFactor}_${i}
  echo "$DIR/create_dataset.R --fit $FIT --data $DATA --subject SubID --covariates $COVARIATES --seed $i --logFC $LOGFC --pDE 0.05 --libScaleFactor ${libScaleFactor} --output $OUTFOLDER/sim_${ID}.h5ad" >> $OUTFOLDER/script_sim.sh
done
done
done

# run sims
cat $OUTFOLDER/script_sim.sh | parallel -P 62

# check that files were written
# cat $OUTFOLDER/script_sim.sh | awk '{print $NF}' | xargs ls > /dev/null 2> file.txt

# cat file.txt | sed -n "s/^ls: cannot access '\(.*\)': No such file or directory$/\1/p" | parallel -P1  basename {} .h5ad | grep -f - $OUTFOLDER/script_sim.sh > jobs.sh

# cat jobs.sh | parallel -P 10


# Concat H5AD
##############

ml purge
ml anaconda3/latest parallel
conda activate jul2026

SRC=/hpc/users/hoffmg01/work/GenomicDataStream_analysis/concat_h5ad.py

echo "" > $OUTFOLDER/script_concat.sh
for N in $(echo $NSAMPLES)
do
for libScaleFactor in $(echo $LSF)
do
for i in $(seq 1 1 $NREPS)
do
  ID=${N}_${libScaleFactor}_${i}
  H5ADS=$(ls sim_${ID}_*.h5ad)
  OUT=sim_${ID}.h5ad
  FILE=file_${ID}.txt
  echo $H5ADS | tr ' ' '\n' > $FILE

  echo "$SRC -i $FILE -o $OUT" >> $OUTFOLDER/script_concat.sh
done
done
done

cat $OUTFOLDER/script_concat.sh | parallel -P 50

# check that files were written
cat $OUTFOLDER/script_concat.sh | awk '{print $NF}' | xargs ls > /dev/null 2> file.txt


# Recode for efficient access
#----------------------------
SRC=/hpc/users/hoffmg01/work/GenomicDataStream_analysis/recode_h5ad.py 

echo "" > $OUTFOLDER/script_recode.sh
for N in $(echo $NSAMPLES)
do
for libScaleFactor in $(echo $LSF)
do
for i in $(seq 1 1 $NREPS)
do
  ID=${N}_${libScaleFactor}_${i}
  FILE=$OUTFOLDER/sim_${ID}.h5ad
  OUT=$OUTFOLDER/sim_${ID}_recode.h5ad
  echo "$SRC --input $FILE --sortBy subclass,SubID --format CSC --compression lzf --out $OUT " >> $OUTFOLDER/script_recode.sh
done
done
done

cat $OUTFOLDER/script_recode.sh | parallel -P 60

# check that files were written
cat $OUTFOLDER/script_recode.sh | awk '{print $NF}' | xargs ls > /dev/null

# remove tmp h5ad's
comm -3 <(ls $OUTFOLDER/*recode.h5ad | sort) <(ls $OUTFOLDER/*.h5ad | sort) | xargs -n1 rm -f

rm -f file*.txt


# Run DE analysis
#----------------

# source ~/.bash_profile

# METHODS=/hpc/users/hoffmg01/work/lucida_analysis/simulations/methods.in
# LOG=$OUTFOLDER/logs
# mkdir -p $LOG

# echo "" > $OUTFOLDER/script_de.sh
# for N in $(echo $NSAMPLES)
# do
# for libScaleFactor in $(echo $LSF)
# do
# for i in $(seq 1 1 $NREPS)
# do
#   ID=${N}_${libScaleFactor}_${i}
#   FILE=$OUTFOLDER/sim_${ID}_recode.h5ad
#   OUT=$OUTFOLDER/res_sim_${ID}.parquet
#   echo "$DIR/run_analysis.R --h5ad $FILE --formula \"~ Dx + (1|SubID)\" --coefTest DxDisease --cluster_id subclass --methods $METHODS --output $OUT 2>&1 > $LOG/${ID}.log " >> $OUTFOLDER/script_de.sh
# done
# done
# done

# cat $OUTFOLDER/script_de.sh | parallel -P 60

# # check that files were written
# cat $OUTFOLDER/script_de.sh | awk '{print $16}' | xargs ls > /dev/null 2> err.log

# cat err.log | tr "'" " " | awk '{print $4}' | parallel basename {} .parquet > jobs.prefix

# grep -f jobs.prefix $OUTFOLDER/script_de.sh | parallel -P 60

# Run Distributed DE
####################

# conda deactivate; source ~/.bash_profile

METHODS=/hpc/users/hoffmg01/work/lucida_analysis/simulations/methods.in
LOG=$OUTFOLDER/logs
mkdir -p $LOG 
mkdir -p $OUTFOLDER/jobs/ 
NTHREADS=12

plugin_dir=$(Rscript -e 'cat(rhdf5filters::hdf5_plugin_path())')

for N in $(echo $NSAMPLES)
do
for libScaleFactor in $(echo $LSF)
do
for i in $(seq 1 1 $NREPS)
do
for METHOD in $(grep -v "#" $METHODS | sed 's/"//g')
do
  ID=${N}_${libScaleFactor}_${i}
  FILE=$OUTFOLDER/sim_${ID}_recode.h5ad
  OUT=$OUTFOLDER/res_sim_${ID}_${METHOD}.parquet
  JOB=$OUTFOLDER/jobs/script_${ID}_${METHOD}.sh

  if [[ ("$METHOD" == "nebula") || ("$METHOD" == "nebula_HL") ]]; then 
    MEM=16000; 
  else 
    MEM=3000; 
  fi

  echo '#!/bin/bash' > $JOB
  echo "#BSUB -P acc_CommonMind
#BSUB -q premium
#BSUB -J ${ID}_${METHOD}
#BSUB -n $NTHREADS
#BSUB -R span[hosts=1] 
#BSUB -R rusage[mem=$MEM]
#BSUB -W 6:00
#BSUB -e $LOG/${ID}_${METHOD}.err
#BSUB -o $LOG/${ID}_${METHOD}.out" >> $JOB
  echo -e "\\nsource ~/.bash_profile\\n" >> $JOB
  
  echo -e "export HDF5_PLUGIN_PATH=$plugin_dir\\n" >> $JOB

  echo "$DIR/run_analysis.R --h5ad $FILE --formula \"~ Dx + (1|SubID) + Age + Sex + PMI\" --coefTest DxDisease --cluster_id subclass --nthreads $NTHREADS --methods $METHOD --output $OUT" >> $JOB
done
done
done
done

# rm -f logs/*
# rm -f res_sim_*


# submit jobs
ls $OUTFOLDER/jobs/* | grep dreamlet | parallel -P1 "bsub < {}"
ls $OUTFOLDER/jobs/* | grep DESeq2 | parallel -P1 "bsub < {}"
ls $OUTFOLDER/jobs/* | grep edgeR | parallel -P1 "bsub < {}"
ls $OUTFOLDER/jobs/* | grep glmGamPoi | parallel -P1 "bsub < {}"

ls $OUTFOLDER/jobs/* | grep lucida | parallel -P1 "bsub < {}"
ls $OUTFOLDER/jobs/* | grep nebula | parallel -P1 "bsub < {}"
ls $OUTFOLDER/jobs/* | grep MAST | parallel -P1 "bsub < {}"

# resub
grep -r AVX512_FP16 logs/*.err | cut -f1 -d':' | cut -f1 -d'.' | xargs -n 1 basename | parallel -P1 echo "jobs/script_{}.sh" | parallel -P1 "bsub < {}"


# ls $OUTFOLDER/jobs/* | grep "dreamlet\|DESeq2\|lucida" | parallel -P1 cat | grep 'analysis' | parallel -P62




# remove logs
# grep -r AVX512_FP16 logs/*.err | cut -f1 -d':' | cut -f1 -d'.' | xargs -n 1 basename | parallel -P1 echo "logs/{}.err" 

ml R/4.6.1


library(Rfast)
X = matrnorm(2000, 10000)
system.time(dcmp <- svd(X))




# Performance plots
###################

ml pandoc
cd /sc/arion/work/hoffmg01/lucida_analysis/simulations

# rm -rf plot_results_PsychAD_cache/ plot_results_PsychAD_files/

rmarkdown::render("plot_results_PsychAD.Rmd")


system("cp -f plot_results_PsychAD.html ~/www/")

https://hoffmg01.dmz.hpc.mssm.edu/plot_results_PsychAD.html



# Memory usage
##############

echo "" > $OUTFOLDER/script_de_mem.sh

for N in $(echo $NSAMPLES)
do
for libScaleFactor in $(echo $LSF)
do
for i in $(seq 1 1 10)
do
for METHOD in $(echo "lucida nebula")
do
  ID=${N}_${libScaleFactor}_${i}
  FILE=$OUTFOLDER/sim_${ID}_recode.h5ad
  OUT=/dev/null
  LOG=$OUTFOLDER/mem_${METHOD}_${ID}.log

  echo "/usr/bin/time -v $DIR/run_analysis.R --h5ad $FILE --formula \"~ Dx + (1|SubID)\" --coefTest DxDisease --cluster_id subclass --methods <(echo \"$METHOD\") --output $OUT 2> $LOG" >> $OUTFOLDER/script_de_mem.sh
done
done
done
done

cat $OUTFOLDER/script_de_mem.sh | parallel -P 35



