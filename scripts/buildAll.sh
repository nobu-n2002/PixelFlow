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
#   This script builds sorce code and creates a new work directory.                 #
# _________________________________________________________________________________ #

echo "\e[1;32m"
echo "           _|                      _|      _|_|  _|                                "
echo " _|_|_|        _|    _|    _|_|    _|    _|      _|    _|_|    _|      _|      _|  "
echo " _|    _|  _|    _|_|    _|_|_|_|  _|  _|_|_|_|  _|  _|    _|  _|      _|      _|  "
echo " _|    _|  _|  _|    _|  _|        _|    _|      _|  _|    _|    _|  _|  _|  _|    "
echo " _|_|_|    _|  _|    _|    _|_|_|  _|    _|      _|    _|_|        _|      _|      "
echo " _|                                                                                "
echo " _|                                                                                "
echo "                                                                                   "
echo " Author: Nobuto NAKAMICHI, Younghwa CHO, Nobuyuki OSHIMA                           "
echo " Date: 25.02.2024                                                                  "
echo " Description:                                                                      "
echo "   This script builds sorce code and creates a new work directory.                 "
echo "\e[0m"

# Function to kill spinner process
cleanup() {
    echo ""
    echo "Termination signal was sent."
    echo "Forcing termination of the build."
    kill $spinner_pid >/dev/null 2>&1
    echo ""
    exit 1
}

# Trap Ctrl+C and clean up
trap cleanup 2

sleep 2

# Start the project build
echo "Starting project build..."

# build 2d program (openMP parallel code for Multi CPUs)
# Start spinner animation in background
(
    while true; do
        printf " / "
        sleep 0.1
        printf "\b\b\b"
        printf " - "
        sleep 0.1
        printf "\b\b\b"
        printf " \\ "
        sleep 0.1
        printf "\b\b\b"
        printf " | "
        sleep 0.1
        printf "\b\b\b"
    done
) &
spinner_pid=$!

# ===============================================================

sh scripts/build2d.sh
sh scripts/build3d.sh

# Stop spinner animation
kill $spinner_pid >/dev/null 2>&1