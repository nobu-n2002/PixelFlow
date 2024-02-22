#!/bin/bash

SRC_DIR=src
SRC=ibm_3d.f90
FC=gfortran
FC_FLAG='-O3 -fopenmp -fno-automatic -o'
BIN_DIR=bin
EXE=ibm3

echo mkdir -p ${BIN_DIR}
echo cd ${BIN_DIR}

mkdir -p ${BIN_DIR}
cd ${BIN_DIR}

echo rm -f *.o *.mod *.exe *.out ibm3
rm -f *.o *.mod *.exe *.out ibm3

echo ${FC} ${FC_FLAG} ${EXE} ../${SRC_DIR}/${SRC}
if ${FC} ${FC_FLAG} ${EXE} ../${SRC_DIR}/${SRC}; then
    echo "Build complete. Executable: bin/ibm3"
else
    echo "Error: Build failed."
fi
