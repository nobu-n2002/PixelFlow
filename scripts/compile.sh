#!/bin/bash
gfortran -O3 -fopenmp -fno-automatic -o bin/ibm3 src/$1
