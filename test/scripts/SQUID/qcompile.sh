#!/bin/bash

module load BaseVEC/2023

SRC_DIR=src
SRC=ibm_3d.f90
#FC=gfortran
#FC=nvfortran
FC=nfort
BIN_DIR=bin
EXE=ibm3

case "$FC" in
#[MEMO] IF "2055 Segmentation fault (core dumped)" APPEARS, ADD "-fno-automatic".
        # gfortran) FC_FLAG='-O3 -fopenmp -fno-automatic -foffload=nvptx-none';;
        gfortran) FC_FLAG='-O3 -fopenmp -fno-automatic -o';;

#[MEMO] IF "Segmentation fault" APPEARS, ADD "-Msave".
        nvfortran) FC_FLAG='-O3 -mp=multicore -Msave -o';;

        nfort) FC_FLAG='-O3 -fopenmp -o'

esac

rm -f ${BIN_DIR}/*.o ${BIN_DIR}/*.mod ${BIN_DIR}/*.exe ${BIN_DIR}/*.out ${BIN_DIR}/ibm3

echo ${FC} ${FC_FLAG} ${BIN_DIR}/${EXE} ${SRC_DIR}/${SRC}

if ${FC} ${FC_FLAG} ${BIN_DIR}/${EXE} ${SRC_DIR}/${SRC}; then
    echo "Build complete. Executable:" ${BIN_DIR}/${EXE}
else
    echo "Error: Build failed."
fi