#!/bin/bash

# Set the DIMENSION variable to 2 for a 2D cylindrical case
# Set the DIMENSION variable to 3 for a 3D Stanford Dragon case
DIMENSION=2

# make directory
mkdir -p bin config etc logs

# target directory
TARGET_DIR="data"

# if target directory is not exist
if [ ! -d "$TARGET_DIR" ]; then
    echo "Directory '$TARGET_DIR' not found. Extracting data.zip..."

    # unzip data.zip
    unzip data.zip

    echo "Extraction complete."
else
    echo "Directory '$TARGET_DIR' already exists. No need to extract."
fi

if [ "$DIMENSION" -eq 3 ]; then
    echo "Making 3d project"
    sh ./scripts/build3d.sh
    sh ./scripts/test3d.sh
elif [ "$DIMENSION" -eq 2 ]; then
    echo "Making 2d project"
    sh ./scripts/build2d.sh
    sh ./scripts/test2d.sh
else
    echo "Invalid value for DIMENSION. Please set DIMENSION to 2 or 3."
    exit 1
fi
