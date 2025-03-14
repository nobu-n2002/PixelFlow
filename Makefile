#########################################################################
#                                                                       #
# PixelFlow :   Flow simulation solver Fluid simulation solver          #
#               implementing the immersed boundary method.              #
#                                                                       #
# Copyright (c) 2023 Nobuto NAKAMICHI, Younghwa CHO, Nobuyuki OSHIMA    #
# All rights reserved.                                                  #
#                                                                       #
#########################################################################

# Define variables
SHELL := /bin/bash
SCRIPTS_DIR := scripts
BUILD_DIR := build 
PROJECT_DIR := projects
f := my_projects

# Function to display usage
usage:
	@echo "Usage: make [OPTIONS]"
	@echo "options:"
	@echo "  build				Build all source codes"
	@echo "  project f=<folder_name>	Create a new folder with the specified name"
	@echo "  help				Display Usage"
	@exit 0

clean:
	@echo "Clean up"
	$(RM) -rf build
	$(RM) -rf bin
	$(RM) -rf lib

# Target for building
build: clean
	cmake -S . -B $(BUILD_DIR)
	cmake --build $(BUILD_DIR)

# Target for creating a new project folder
project:
	@echo ""
	@echo "------------------------------------------"
	@echo " Making project $(f)"
	@echo "------------------------------------------"
	@echo ""
	@if [ ! -d "$(PROJECT_DIR)/$(f)" ]; then \
	    echo "Creating folder $(f)"; \
	    mkdir -p "$(PROJECT_DIR)/$(f)"; \
	    cp -r template/* "$(PROJECT_DIR)/$(f)"; \
	    echo ""; \
	    echo "Project $(f) was completed to build."; \
	else \
	    echo "Error: $(PROJECT_DIR)/$(f) already exists."; \
	fi

# Target for displaying usage
help:
	$(MAKE) usage

testproject:
	$(MAKE) build
	$(MAKE) project f=cylinder-2d
	${RM} $(PROJECT_DIR)/cylinder-2d/config/controlDict.txt
	${RM} $(PROJECT_DIR)/cylinder-2d/config/omp_config
	@unzip test/cylinder-2d.zip -d $(PROJECT_DIR)
	$(MAKE) project f=stanford-dragon
	${RM} $(PROJECT_DIR)/stanford-dragon/config/controlDict.txt
	${RM} $(PROJECT_DIR)/stanford-dragon/config/omp_config
	@unzip test/stanford-dragon.zip -d $(PROJECT_DIR)
	$(MAKE) project f=room
	${RM} $(PROJECT_DIR)/room/config/controlDict.txt
	${RM} $(PROJECT_DIR)/room/config/omp_config
	@unzip test/room.zip -d $(PROJECT_DIR)
