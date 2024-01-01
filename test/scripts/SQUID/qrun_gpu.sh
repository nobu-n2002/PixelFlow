#!/bin/bash
#------- qsub option -----------
#PBS -q SQUID
#PBS --group=hp000000                           # Group ID
#PBS -l elapstim_req=05:00:00                   # Forcibly terminating after 5 hours.
#PBS -l gpunum_job=1
#PBS -m eb                                      # At the start of the batch request execution, send an email.
#PBS -M hoge@mail                               # e-mail address
#------- Program execution -----------
module purge
module load BaseGPU
export UCX_TLS=sm,cuda_copy,cuda_ipc,gdr_copy,self  #1ノードで実行する場合は設定してください

cd $PBS_O_WORKDIR

EXE=ibm2
LOG_DIR=logs

STDOUT_FNAME=runlog_$(date "+%Y.%m.%d-%H.%M.%S").txt 
echo 'GPU used for this run '> $LOG_DIR/$STDOUT_FNAME
./bin/$EXE >> $LOG_DIR/$STDOUT_FNAME 2>&1 &