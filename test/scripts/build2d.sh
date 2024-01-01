#!/bin/bash

SRC_DIR=src
SRC=ibm_2d.f90
FC=gfortran
FC_FLAG='-O3 -fopenmp -fno-automatic -o'
BIN_DIR=bin
EXE=ibm2

mkdir -p bin

rm -f ${BIN_DIR}/*.o ${BIN_DIR}/*.mod ${BIN_DIR}/*.exe ${BIN_DIR}/*.out ${BIN_DIR}/ibm2

echo ${FC} ${FC_FLAG} ${BIN_DIR}/${EXE} ${SRC_DIR}/${SRC}

if ${FC} ${FC_FLAG} ${BIN_DIR}/${EXE} ${SRC_DIR}/${SRC}; then
    echo "Build complete. Executable: bin/ibm2"
else
    echo "Error: Build failed."
fi