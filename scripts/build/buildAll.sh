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
#   This script builds sorce code and creates a new work directory.                 #
# _________________________________________________________________________________ #

DIR=scripts/build

sh ${DIR}/clean.sh
sh ${DIR}/build_omp.sh -f ibm_2d_uniform_omp_cpu.f90 -o ibm2_uniform_omp
sh ${DIR}/build_omp.sh -f ibm_2d_drag_omp_cpu.f90 -o ibm2_drag_omp
sh ${DIR}/build_omp.sh -f ibm_2d_backstep_omp_cpu.f90 -o ibm2_backstep_omp
sh ${DIR}/build_omp.sh -f ibm_3d_uniform_omp_cpu.f90 -o ibm3_uniform_omp
sh ${DIR}/build_omp.sh -f ibm_3d_air_condition_omp_cpu.f90 -o ibm3_air_condition_omp
sh ${DIR}/build_acc.sh -f ibm_2d_uniform_acc_gpu.f90 -o ibm2_uniform_acc
sh ${DIR}/build_acc.sh -f ibm_2d_drag_acc_gpu.f90 -o ibm2_drag_acc
sh ${DIR}/build_acc.sh -f ibm_3d_uniform_acc_gpu.f90 -o ibm3_uniform_acc
sh ${DIR}/build_acc.sh -f ibm_3d_air_condition_acc_gpu.f90 -o ibm3_air_condition_acc
