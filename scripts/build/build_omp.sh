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
BIN_DIR=bin

# echo mkdir -p ${BIN_DIR}
# echo cd ${BIN_DIR}
mkdir -p ${BIN_DIR}
cd ${BIN_DIR}

# default (dummy)
SRC=main.f90
EXE=main

# compiler flag
FC=gfortran
FC_FLAG='-O3 -fopenmp -fno-automatic'

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
    echo ${FC} ${FC_FLAG} ../${SRC_DIR}/${SRC} -o ${EXE} 
    if ${FC} ${FC_FLAG} ../${SRC_DIR}/${SRC} -o ${EXE}; then
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