#!/bin/bash
#------- qsub option -----------
#PBS -q SQUID
#PBS --group=hp000000                           # Group ID
#PBS -l elapstim_req=05:00:00                   # Forcibly terminating after 5 hours.
#PBS --venode=1                                 # Using one node for vector computation.
#PBS -m eb                                      # At the start of the batch request execution, send an email.
#PBS -M hoge@mail                               # e-mail address
#------- Program execution -----------
module purge
module load BaseVEC/2023

cd $PBS_O_WORKDIR

STDOUT_FNAME=runlog_$(date "+%Y.%m.%d-%H.%M.%S").txt
echo 'Use VECTOR machine' > $STDOUT_FNAME
echo 'Number of threads used = '$OMP_NUM_THREADS > $STDOUT_FNAME
./bin/ibm3 controlDict.txt >>$STDOUT_FNAME 2>&1 &