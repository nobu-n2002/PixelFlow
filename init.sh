#!/bin/bash

# Make directory
mkdir -p bin data config etc

# Build src *.f90 and create controlDict.txt
echo "Making 3d project"
sh ./scripts/build3d.sh
echo "Making 2d project"
sh ./scripts/build2d.sh
sh ./scripts/project.sh
