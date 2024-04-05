# Enable verbose output
# set(CMAKE_VERBOSE_MAKEFILE ON)

# Set directories
set(BIN_DIR ${CMAKE_CURRENT_SOURCE_DIR}/../../bin)
set(SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR})
set(LIB_DIR ${SRC_DIR}/../../lib/modules)

# Compile library files
file(STRINGS ${CMAKE_CURRENT_SOURCE_DIR}/mainFileLists.txt MAIN_FILES)

# Compile each source file separately
foreach(SRC ${MAIN_FILES})

    # Get executable file name
    get_filename_component(EXE_NAME ${SRC} NAME_WE)

    # Create executable
    add_executable(${EXE_NAME} ${SRC})

    # Set output directory
    set_target_properties(${EXE_NAME} PROPERTIES
        RUNTIME_OUTPUT_DIRECTORY ${BIN_DIR}
    )

    # Include directories
    target_include_directories(${EXE_NAME} PRIVATE ${LIB_DIR})

    # Link libraries
    target_link_libraries(${EXE_NAME}
        PRIVATE
        ibm_library
    )

endforeach()