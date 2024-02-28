#!/bin/bash

# Function to display script usage
display_usage() {
    echo "Usage: $0 [OPTIONS] <folder_name>"
    echo "Options:"
    echo "  -b                 Build all source codes"
    echo "  -f <folder_name>   Create a new folder with the specified name"
    echo "  -h                 Display Usage"
}

# Initialize flags
run_build=false
create_folder=false
folder_name=""

# Parse options
while [ "$#" -gt 0 ]; do
    case $1 in
        -b) run_build=true;;
        -build) run_build=true;;
        -f) create_folder=true; folder_name="$2"; shift ;;
        -folder) create_folder=true; folder_name="$2"; shift ;;
        -bf|-fb) run_build=true; create_folder=true; folder_name="$2"; shift ;;
        -h|--help) display_usage; exit 0 ;;
        *) echo "Error: Invalid option '$1'"; display_usage; exit 1 ;;
    esac
    shift
done

# Check if both -b and -f options are provided
if [ "$run_build" = true ]; then
    # Run buildAll.sh script
    sh scripts/buildAll.sh
fi

# Check if -f option is provided
if [ "$create_folder" = true ]; then
    # Create the folder
    mkdir -p "$folder_name"

    # Check if folder creation was successful
    if [ $? -eq 0 ]; then
        echo "Folder '$folder_name' created successfully."
    else
        echo "Error: Failed to create folder '$folder_name'."
        exit 1
    fi

    # Make directory
    cd "$folder_name" || exit
    echo "Making directory"
    echo "mkdir -p data config"
    mkdir -p data config

    echo "../scripts/control.sh"
    sh ../scripts/control.sh

    echo ""
    echo "Project $folder_name was completed to build."
    echo ""
fi