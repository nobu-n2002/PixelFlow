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
#   This script builds 2d incompressible flow simulation code.                      #
# _________________________________________________________________________________ #

# Function to kill spinner process
cleanup() {
    echo ""
    echo "Termination signal was sent."
    echo "Forcing termination of the build."
    kill $spinner_pid >/dev/null 2>&1
    echo ""
    exit 1
}

# Trap Ctrl+C and clean up
trap cleanup 2

sleep 2

# Start the project build
echo ""
echo "Making 2d project..."

SRC_DIR=src
BIN_DIR=bin

# echo mkdir -p ${BIN_DIR}
# echo cd ${BIN_DIR}
mkdir -p ${BIN_DIR}
cd ${BIN_DIR}

# build 2d program (openMP parallel code for Multi CPUs)
# Start spinner animation in background
(
    while true; do
        printf " / "
        sleep 0.1
        printf "\b\b\b"
        printf " - "
        sleep 0.1
        printf "\b\b\b"
        printf " \\ "
        sleep 0.1
        printf "\b\b\b"
        printf " | "
        sleep 0.1
        printf "\b\b\b"
    done
) &
spinner_pid=$!

# ===============================================================
echo ""
echo 'Build 2d program (openMP parallel code for Multi CPUs)'
SRC1=ibm_2d_omp_cpu.f90
SRC2=ibm_2d_drag_omp_cpu.f90
FC=gfortran
FC_FLAG='-O3 -fopenmp -fno-automatic -o'
EXE1=ibm2_omp
EXE2=ibm2_drag_omp

echo rm -f *.o *.mod *.exe *.out ${EXE1} ${EXE2}
rm -f *.o *.mod *.exe *.out ${EXE1}  ${EXE2}

# Check if GNU compiler is available
if [ -x "$(command -v gfortran)" ]; then
    echo "--- gfortran is installed. Compiling with gfortran..."
    echo ${FC} ${FC_FLAG} ${EXE1} ../${SRC_DIR}/${SRC1}
    if ${FC} ${FC_FLAG} ${EXE1} ../${SRC_DIR}/${SRC1}; then
        echo "--- Build complete. Executable: bin/${EXE1}"
    else
        echo "\e[1;31mError:\e[0m Build failed."
    fi
    echo ""
    echo ${FC} ${FC_FLAG} ${EXE2} ../${SRC_DIR}/${SRC2}
    if ${FC} ${FC_FLAG} ${EXE2} ../${SRC_DIR}/${SRC2}; then
        echo "--- Build complete. Executable: bin/${EXE2}"
    else
        echo "\e[1;31mError:\e[0m Build failed."
    fi
else
    echo "\e[1;35mWarning:\e[0m gfortran is not installed. If you want to use GPU-accelerated compilation, please download nvfortran."
    
fi

# Stop spinner animation
kill $spinner_pid >/dev/null 2>&1

# ===============================================================
# build 2d program (openACC parallel code for GPU device)
echo ""
echo 'Build 2d program (openACC parallel code for GPU device)'
SRC1=ibm_2d_acc_gpu.f90
SRC2=ibm_2d_drag_acc_gpu.f90
FC=nvfortran
FC_FLAG='-O3 -acc=gpu -gpu=ccall -Minfo=accel -o'
EXE1=ibm2_acc
EXE2=ibm2_drag_acc

echo rm -f *.o *.mod *.exe *.out ${EXE1} ${EXE2}
rm -f *.o *.mod *.exe *.out ${EXE1}  ${EXE2}

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
            echo "\e[1;31mError:\e[0m Build failed."
        fi
        echo ${FC} ${FC_FLAG} ${EXE2} ../${SRC_DIR}/${SRC2}
        if ${FC} ${FC_FLAG} ${EXE2} ../${SRC_DIR}/${SRC2}; then
            echo "Build complete. Executable: bin/${EXE2}"
        else
            echo "\e[1;31mError:\e[0m Build failed."
        fi
    else
        echo "\e[1;35mWarning:\e[0m nvfortran is not installed. If you want to use GPU-accelerated compilation, please download nvfortran."
    fi
else
    echo "\e[1;35mWarning:\e[0m No GPU devices available. Unable to compile GPU-accelerated code."
    
fi

# Stop spinner animation
kill $spinner_pid >/dev/null 2>&1

# go to PixelFlow/
cd ..

echo ""
echo "2d program build complete."
echo ""

exit 0