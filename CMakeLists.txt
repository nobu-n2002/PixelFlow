#########################################################################
#                                                                       #
# PixelFlow :   Flow simulation solver Fluid simulation solver          #
#               implementing the immersed boundary method.              #
#                                                                       #
# Copyright (c) 2023 Nobuto NAKAMICHI, Younghwa CHO, Nobuyuki OSHIMA    #
# All rights reserved.                                                  #
#                                                                       #
#########################################################################

cmake_minimum_required(VERSION 3.12)
project(IBM
    VERSION 1.0
    DESCRIPTION "Flow Simulation Code implementing IB Method"
    LANGUAGES Fortran C CXX)

# Build option setting
option(with_openmp "Enable OpenMP" ON)
option(with_openacc "Enable OpenACC" ON)

set(DCMAKE_BUILD_TYPE "Release")

if(with_openmp)
    find_package(OpenMP REQUIRED)
    if(OpenMP_FOUND)
        message(STATUS "OpenMP found")
        # Build libraries
        add_subdirectory(src/omp_parallel/lib)
        # Build source file with omp
        add_subdirectory(src/omp_parallel)
        # add_subdirectory(src/dev)
    else()
        message(FATAL_ERROR "OpenMP not found")
    endif()
endif()

if(with_openacc)
    find_package(OpenACC REQUIRED)
    # Check for CUDA and nvfortran
    find_package(CUDA QUIET)
    if(CUDA_FOUND)
        message(STATUS "CUDA found.")
        # Check if nvfortran is available
        include(CheckLanguage)
        check_language(CUDA)
        if(CMAKE_CUDA_COMPILER)
            message(STATUS "nvfortran found.")
            # Build libraries
            add_subdirectory(src/oacc_parallel/lib)
            # Build source file with oacc
            add_subdirectory(src/oacc_parallel)
        else()
            message(WARNING "nvfortran not found.")
        endif()
        
    else()
        message(WARNING "CUDA not found. GPU acceleration may not be available.")
    endif()
endif()