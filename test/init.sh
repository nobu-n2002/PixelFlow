#!/bin/bash

DIMENSION=3

# Make directory
mkdir -p bin data config etc

unzip data.zip

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
