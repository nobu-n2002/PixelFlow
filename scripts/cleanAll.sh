#!/bin/bash

# _________________________________________________________________________________ #
#                                                                                   #
#           _|                      _|      _|_|  _|                                #
# _|_|_|        _|    _|    _|_|    _|    _|      _|    _|_|    _|      _|      _|  #
# _|    _|  _|    _|_|    _|_|_|_|  _|  _|_|_|_|  _|  _|    _|  _|      _|      _|  #
# _|    _|  _|  _|    _|  _|        _|    _|      _|  _|    _|    _|  _|  _|  _|    #
# _|_|_|    _|  _|    _|    _|_|_|  _|    _|      _|    _|_|        _|      _|      #
# _|                                                                                #
# _|                                                                                #
#                                                                                   #
# Author: Nobuto NAKAMICHI, Younghwa CHO, Nobuyuki OSHIMA                           #
# Date: 25.02.2024                                                                  #
# Description:                                                                      #
#   This script clean all runlogs and output files without output_folder/*.vtk      #
# _________________________________________________________________________________ #

while true; do
    echo "Are you sure you want to delete all output files? (yes/no)"
    read answer

    case "$answer" in
        [yY]|[yY][eE][sS])
            echo "Deleting files..."
            rm -f logs/runlog* etc/*
            echo "Files have been deleted."
            break
            ;;
        [nN]|[nN][oO])
            echo "Deletion canceled."
            break
            ;;
        *)
            echo "Invalid input. Please enter 'yes' or 'no'."
            ;;
    esac
done