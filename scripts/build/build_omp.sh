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
#   This script builds OpenMP parallel code for Multi CPUs.                      #
# _________________________________________________________________________________ #

SRC_DIR=src
SRC_OMP_DIR=${SRC_DIR}/omp_parallel
BIN_DIR=bin
LIB_DIR=${SRC_DIR}/lib

# echo mkdir -p ${BIN_DIR}
# echo cd ${BIN_DIR}
mkdir -p ${BIN_DIR}
cd ${BIN_DIR}

# default (dummy)
SRC=omp_parallel/main.f90
EXE=main

# Library and Shared Files
LIB1=global.f90
LIB2=utils.f90
LIB3=output.f90

# compiler flag
FC=gfortran
FC_FLAG='-O3 -fopenmp -fno-automatic -mcmodel=medium'
# FC_FLAG_DBG='-Wall'

while getopts ":f:o:" opt; do
  case ${opt} in
    f )
      SRC=$OPTARG
      ;;
    o )
      EXE=$OPTARG
      ;;
    \? )
      echo "Invalid option: $OPTARG" 1>&2
      exit 1
      ;;
    : )
      echo "Option -$OPTARG requires an argument." 1>&2
      exit 1
      ;;
  esac
done
shift $((OPTIND -1))

echo ""
echo '------------------------------------------'
echo 'Build OpenMP parallel code for Multi CPUs'
echo 'source file:' ${SRC}
echo 'executable file:' ${EXE}
echo '------------------------------------------'

echo rm -f *.o *.mod *.out ${EXE} ${EXE}.exe
rm -f *.o *.mod *.out ${EXE} ${EXE}.exe

# Check if GNU compiler is available
if [ -x "$(command -v gfortran)" ]; then
    echo "--- gfortran is installed. Compiling with gfortran..."

    # Compiles but does not link
    echo ${FC} ${FC_FLAG} ${FC_FLAG_DBG} -c ../${LIB_DIR}/${LIB1}
    ${FC} ${FC_FLAG} ${FC_FLAG_DBG} -c ../${LIB_DIR}/${LIB1}
    
    echo ${FC} ${FC_FLAG} ${FC_FLAG_DBG} -c ../${LIB_DIR}/${LIB2}
    ${FC} ${FC_FLAG} ${FC_FLAG_DBG} -c ../${LIB_DIR}/${LIB2}
    
    echo ${FC} ${FC_FLAG} ${FC_FLAG_DBG} -c ../${LIB_DIR}/${LIB3}
    ${FC} ${FC_FLAG} ${FC_FLAG_DBG} -c ../${LIB_DIR}/${LIB3}

    echo ${FC} ${FC_FLAG} ${FC_FLAG_DBG} -c ../${SRC_OMP_DIR}/${SRC}
    ${FC} ${FC_FLAG} ${FC_FLAG_DBG} -c ../${SRC_OMP_DIR}/${SRC}
    
    # Link all *.o files in the bin folder
    echo ${FC} ${FC_FLAG} ${FC_FLAG_DBG} -o ${EXE} *.o
    if ${FC} ${FC_FLAG} ${FC_FLAG_DBG} *.o -o ${EXE} ; then
        echo "--- Build complete. Executable: bin/${EXE}"
    else
        printf "\e[1;31mError:\e[0m Build failed.\\n"
    fi
else
    printf "\e[1;35mWarning:\e[0m gfortran is not installed. If you want to use GNU-compilation, please download gfortran.\\n"
    printf "\e[1;31mError:\e[0m Build failed.\\n"
fi

# go back to PixelFlow/
cd ..

exit 0