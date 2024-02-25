#!/bin/bash

# Check if an argument is provided
if [ $# -ne 1 ]; then
    echo "Usage: $0 <folder_name>"
    exit 1
fi

# Get folder name from the argument
folder_name="$1"

# Create the folder
mkdir "$folder_name"

# Check if folder creation was successful
if [ $? -eq 0 ]; then
    echo "Folder '$folder_name' created successfully."
else
    echo "\e[1;31m Error: \e[0m Failed to create folder '$folder_name'."
    exit 1
fi

sh scripts/buildAll.sh

# Make directory
cd "$folder_name"
echo "Making directory"
echo "mkdir -p data config"
mkdir -p data config

echo "../scripts/control.sh"
sh ../scripts/control.sh

echo ""
echo "Project $folder_name was completed to build."
echo ""
exit 0
