# Enable verbose output
# set(CMAKE_VERBOSE_MAKEFILE ON)

# Set directories
set(SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR})
set(LIB_DIR_NVF ${SRC_DIR}/../../../lib/nvf)
set(CMAKE_Fortran_COMPILER "nvfortran")

# Compile library files
file(GLOB LIB_SOURCES *.f90)

# Set output directory for module files
set(MODULE_OUTPUT_NVF_DIR ${LIB_DIR_NVF}/modules)

# Create library
add_library(ibm_library_nvf ${LIB_SOURCES})

# Set output directories
set_target_properties(ibm_library_nvf
    PROPERTIES
    ARCHIVE_OUTPUT_DIRECTORY ${LIB_DIR_NVF}
    # Fortran_MODULE_DIRECTORY ${MODULE_OUTPUT_NVF_DIR} # Specify directory for module files
)

# Add compiler flags for OpenACC support
target_compile_options(ibm_library_nvf PRIVATE
    $<$<COMPILE_LANGUAGE:Fortran>:-acc=gpu -gpu=ccall -Minfo=accel -mcmodel=medium>
)

# Custom command to copy module files after build
add_custom_command(TARGET ibm_library_nvf POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E make_directory ${MODULE_OUTPUT_NVF_DIR}
    COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_BINARY_DIR}/*.mod ${MODULE_OUTPUT_NVF_DIR}
)
