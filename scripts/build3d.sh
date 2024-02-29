#!/bin/bash

# _________________________________________________________________________________ #
#                                                                                   #
#           _|                      _|      _|_|  _|                                #
# _|_|_|        _|    _|    _|_|    _|    _|      _|    _|_|    _|      _|      _|  #
# _|    _|  _|    _|_|    _|_|_|_|  _|  _|_|_|_|  _|  _|    _|  _|      _|      _|  #
# _|    _|  _|  _|    _|  _|        _|    _|      _|  _|    _|    _|  _|  _|  _|    #
# _|_|_|    _|  _|    _|    _|_|_|  _|    _|      _|    _|_|        _|      _|      #
# _|                                                                                #
# _|                                                                                #
#                                                                                   #
# Author: Nobuto NAKAMICHI, Younghwa CHO, Nobuyuki OSHIMA                           #
# Date: 25.02.2024                                                                  #
# Description:                                                                      #
#   This script builds 3d incompressible flow simulation code.                      #
# _________________________________________________________________________________ #

# Start the project build
echo ""
echo "Making 3d project..."

SRC_DIR=src
BIN_DIR=bin

# echo mkdir -p ${BIN_DIR}
# echo cd ${BIN_DIR}
mkdir -p ${BIN_DIR}
cd ${BIN_DIR}
echo

# ===============================================================
# build 2d program (openMP parallel code for Multi CPUs)
echo ""
echo 'Build 3d program (openMP parallel code for Multi CPUs)'
SRC1=ibm_3d_omp_cpu.f90
SRC2=ibm_3d_air_condition_omp_cpu.f90
FC=gfortran
FC_FLAG='-O3 -fopenmp -fno-automatic -o'
EXE1=ibm3_omp
EXE2=ibm3_air_condition_omp

echo rm -f *.o *.mod *.out ${EXE1} ${EXE2} ${EXE1}.exe ${EXE2}.exe
rm -f *.o *.mod *.out ${EXE1} ${EXE2} ${EXE1}.exe ${EXE2}.exe

# Check if GNU compiler is available
if [ -x "$(command -v gfortran)" ]; then
    echo "--- gfortran is installed. Compiling with gfortran..."
    echo ${FC} ${FC_FLAG} ${EXE1} ../${SRC_DIR}/${SRC1}
    if ${FC} ${FC_FLAG} ${EXE1} ../${SRC_DIR}/${SRC1}; then
        echo "--- Build complete. Executable: bin/${EXE1}"
    else
        printf "\e[1;31mError:\e[0m Build failed.\\n"
    fi
    echo ""
    echo ${FC} ${FC_FLAG} ${EXE2} ../${SRC_DIR}/${SRC2}
    if ${FC} ${FC_FLAG} ${EXE2} ../${SRC_DIR}/${SRC2}; then
        echo "--- Build complete. Executable: bin/${EXE2}"
    else
        printf "\e[1;31mError:\e[0m Build failed.\\n"
    fi
else
    printf "\e[1;35mWarning:\e[0m gfortran is not installed. If you want to use GPU-accelerated compilation, please download nvfortran.\\n"
fi

# ===============================================================
# build 3d program (openACC parallel code for GPU device)
echo ""
echo 'Build 3d program (openACC parallel code for GPU device)'
SRC1=ibm_3d_acc_gpu.f90
SRC2=ibm_3d_air_condition_acc_gpu.f90
FC=nvfortran
FC_FLAG='-O3 -acc=gpu -gpu=ccall -Minfo=accel -o'
EXE1=ibm3_acc
EXE2=ibm3_air_condition_acc

echo rm -f *.o *.mod *.out ${EXE1} ${EXE2} ${EXE1}.exe ${EXE2}.exe
rm -f *.o *.mod *.out ${EXE1} ${EXE2} ${EXE1}.exe ${EXE2}.exe

# Check if GPU devices are available
if [ -x "$(command -v nvidia-smi)" ]; then
    echo "GPU devices are available."
    
    # Check if nvfortran is installed
    if [ -x "$(command -v nvfortran)" ]; then
        echo "nvfortran is installed. Compiling with nvfortran..."
        echo ${FC} ${FC_FLAG} ${EXE1} ../${SRC_DIR}/${SRC1}
        if ${FC} ${FC_FLAG} ${EXE1} ../${SRC_DIR}/${SRC1}; then
            echo "Build complete. Executable: bin/${EXE1}"
        else
            printf "\e[1;31mError:\e[0m Build failed.\\n"
        fi
        echo ${FC} ${FC_FLAG} ${EXE2} ../${SRC_DIR}/${SRC2}
        if ${FC} ${FC_FLAG} ${EXE2} ../${SRC_DIR}/${SRC2}; then
            echo "Build complete. Executable: bin/${EXE2}"
        else
            printf "\e[1;31mError:\e[0m Build failed.\\n"
        fi
    else
        printf "\e[1;35mWarning:\e[0m nvfortran is not installed. If you want to use GPU-accelerated compilation, please download nvfortran.\\n"
    fi
else
    printf "\e[1;35mWarning:\e[0m No GPU devices available. Unable to compile GPU-accelerated code.\\n"
fi

# go to PixelFlow/
cd ..

echo ""
echo "3d program build complete."
echo ""

exit 0
