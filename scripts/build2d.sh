#!/bin/bash

SRC_DIR=src
SRC=ibm_2d.f90
FC=gfortran
FC_FLAG='-O3 -fopenmp -fno-automatic -o'
BIN_DIR=bin
EXE=ibm2

mkdir -p bin
cd bin

rm -f *.o *.mod *.exe *.out ibm2

echo ${FC} ${FC_FLAG} ${BIN_DIR}/${EXE} ${SRC_DIR}/${SRC}

if ${FC} ${FC_FLAG} ${EXE} ../${SRC_DIR}/${SRC}; then
    echo "Build complete. Executable: bin/ibm2"
else
    echo "Error: Build failed."
fi
