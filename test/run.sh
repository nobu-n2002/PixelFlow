#!/bin/bash

# Set the DIMENSION variable to 2 for a 2D cylindrical case
# Set the DIMENSION variable to 3 for a 3D Stanford Dragon case
DIMENSION=2

export OMP_NUM_THREADS=8

STDOUT_FNAME="runlog_$(date "+%Y.%m.%d-%H.%M.%S").txt"
mkdir -p logs
echo "Number of threads used = $OMP_NUM_THREADS" > "logs/$STDOUT_FNAME"

# DIMENSIONの値によってibm3またはibm2を選択して実行
if [ "$DIMENSION" -eq 3 ]; then
    echo "Running ibm3..."
    ./bin/ibm3 >> "logs/$STDOUT_FNAME" 2>&1 &
elif [ "$DIMENSION" -eq 2 ]; then
    echo "Running ibm2..."
    ./bin/ibm2 >> "logs/$STDOUT_FNAME" 2>&1 &
else
    echo "Invalid value for DIMENSION. Please set DIMENSION to 2 or 3."
    exit 1
fi