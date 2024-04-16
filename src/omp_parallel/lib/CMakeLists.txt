# Enable verbose output
# set(CMAKE_VERBOSE_MAKEFILE ON)

# Set directories
set(SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR})
set(LIB_DIR ${SRC_DIR}/../../../lib/gf)

# Compile library files
file(GLOB LIB_SOURCES *.f90)

# Set output directory for module files
set(MODULE_OUTPUT_DIR ${LIB_DIR}/modules)

# Create library
add_library(ibm_library ${LIB_SOURCES})

# Set output directories
set_target_properties(ibm_library
    PROPERTIES
    ARCHIVE_OUTPUT_DIRECTORY ${LIB_DIR}
    Fortran_MODULE_DIRECTORY ${MODULE_OUTPUT_DIR} # Specify directory for module files
)