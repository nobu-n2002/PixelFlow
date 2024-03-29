cmake_minimum_required(VERSION 3.12)
project(STLToPorosityConverter
    VERSION 1.0
    DESCRIPTION "Flow Smulation Code implementing IB Method"
    LANGUAGES CXX)
add_subdirectory(src/omp_parallel)
add_subdirectory(src/oacc_parallel)