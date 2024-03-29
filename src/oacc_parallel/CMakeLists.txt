cmake_minimum_required(VERSION 3.12)
project(OpenACC_Parallel_Code
    VERSION 1.0
    DESCRIPTION "Builds OpenACC parallel code for GPU"
    LANGUAGES Fortran)

# Enable verbose output
# set(CMAKE_VERBOSE_MAKEFILE ON)

# For OpenACC
find_package(OpenACC REQUIRED)

# Check for CUDA and nvfortran
find_package(CUDA)
if(CUDA_FOUND)
    message(STATUS "CUDA found.")
    
    # Check if nvfortran is available
    include(CheckLanguage)
    check_language(CUDA)
    if(CMAKE_CUDA_COMPILER)
        message(STATUS "nvfortran found.")
        set(flag  1)
    else()
        message(WARNING "nvfortran not found.")
        set(flag  0)
    endif()
    
else()
    message(WARNING "CUDA not found. GPU acceleration may not be available.")
    set(flag  0)
endif()

if(flag)
# Set directories
set(BIN_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../../bin)
set(SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR})
set(LIB_DIR ${SRC_DIR}/../lib)

# Compile library files
file(STRINGS ${CMAKE_CURRENT_SOURCE_DIR}/mainFileLists.txt MAIN_FILES)

# Compile each source file separately
foreach(SRC ${MAIN_FILES})
    get_filename_component(EXE_NAME ${SRC} NAME_WE)

    # Create executable
    add_executable(${EXE_NAME} ${SRC})

    # Set output directory
    set_target_properties(${EXE_NAME} PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY ${BIN_DIR}
        COMPILE_FLAGS "${OpenACC_Fortran_FLAGS} ${OpenACC_C_FLAGS}"
        LINK_FLAGS "${OpenACC_Fortran_FLAGS} ${OpenACC_C_FLAGS}"
    )

    # Include directories
    target_include_directories(${EXE_NAME} PRIVATE ${LIB_DIR})

    # Compile library files for each executable
    file(STRINGS ${LIB_DIR}/libFileLists.txt LIB_FILES)

    # Prepend LIB_DIR to each file name in LIB_FILES list
    foreach(LIB_FILE ${LIB_FILES})
        list(APPEND LIB_FILES_WITH_PATH "${LIB_DIR}/${LIB_FILE}")
    endforeach()

    # Use modified LIB_FILES_WITH_PATH as it contains full file paths
    add_library(${EXE_NAME}_LIBRARY_FILES OBJECT ${LIB_FILES_WITH_PATH})

    # Set compiler flags
    set_target_properties(${EXE_NAME}_LIBRARY_FILES PROPERTIES
        COMPILE_FLAGS "${OpenACC_Fortran_FLAGS}"
        LINK_FLAGS "${OpenACC_Fortran_FLAGS}"
    )

    # Link libraries
    target_link_libraries(${EXE_NAME}
        PRIVATE
        ${EXE_NAME}_LIBRARY_FILES
        OpenACC::OpenACC_Fortran
    )
endforeach()
endif()
