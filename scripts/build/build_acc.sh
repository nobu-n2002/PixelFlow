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
#   This script builds  openACC parallel code for GPU device.                      #
# _________________________________________________________________________________ #

SRC_DIR=src
BIN_DIR=bin

# echo mkdir -p ${BIN_DIR}
# echo cd ${BIN_DIR}
mkdir -p ${BIN_DIR}
cd ${BIN_DIR}

# build openACC parallel code for GPU device

# default (dummy)
SRC=main.f90
EXE=main

# compiler flag
FC=nvfortran
FC_FLAG='-O3 -acc=gpu -gpu=ccall -Minfo=accel'

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
echo 'Build openACC parallel code for GPU device'
echo 'source file:' ${SRC}
echo 'executable file:' ${EXE}
echo '------------------------------------------'

echo rm -f *.o *.mod *.out ${EXE} ${EXE}.exe
rm -f *.o *.mod *.out ${EXE} ${EXE}.exe

# Check if GPU devices are available
if [ -x "$(command -v nvidia-smi)" ]; then
    echo "GPU devices are available."
    
    # Check if nvfortran is installed
    if [ -x "$(command -v nvfortran)" ]; then
        echo "nvfortran is installed. Compiling with nvfortran..."
        echo ${FC} ${FC_FLAG} -o ${EXE} ../${SRC_DIR}/${SRC}
        if ${FC} ${FC_FLAG} -o ${EXE} ../${SRC_DIR}/${SRC}; then
            echo "Build complete. Executable: bin/${EXE}"
        else
            printf "\e[1;31mError:\e[0m Build failed.\\n"
        fi
    else
        printf "\e[1;35mWarning:\e[0m nvfortran is not installed. If you want to use GPU-accelerated compilation, please download nvfortran.\\n"
    fi
else
    printf "\e[1;35mWarning:\e[0m No GPU devices available. Unable to compile GPU-accelerated code.\\n"
    printf "\e[1;31mError:\e[0m Build failed.\\n"
fi

# go back to PixelFlow/
cd ..

exit 0