#!/bin/bash
#------- qsub option -----------
#PBS -q SQUID
#PBS --group=hp000000                           # Group ID
#PBS -l elapstim_req=05:00:00                   # Forcibly terminating after 5 hours.
#PBS -l cpunum_job=76
#PBS -v OMP_NUM_THREADS=76
#PBS -m eb                                      # バッチリクエスト実行開始時にメールを送信
#PBS -M hoge@mail                               # e-mail address
#------- Program execution -----------
module purge
module load BaseCPU/2023

cd $PBS_O_WORKDIR

EXE=ibm3
LOG_DIR=logs

STDOUT_FNAME=runlog_$(date "+%Y.%m.%d-%H.%M.%S").txt
echo 'Number of threads used = '$OMP_NUM_THREADS > $LOG_DIR/$STDOUT_FNAME
./bin/$EXE >> $LOG_DIR/$STDOUT_FNAME 2>&1 &