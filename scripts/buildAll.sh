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

printf "\e[1;32m                                                                           \\n"
printf "           _|                      _|      _|_|  _|                                \\n"
printf " _|_|_|        _|    _|    _|_|    _|    _|      _|    _|_|    _|      _|      _|  \\n"
printf " _|    _|  _|    _|_|    _|_|_|_|  _|  _|_|_|_|  _|  _|    _|  _|      _|      _|  \\n"
printf " _|    _|  _|  _|    _|  _|        _|    _|      _|  _|    _|    _|  _|  _|  _|    \\n"
printf " _|_|_|    _|  _|    _|    _|_|_|  _|    _|      _|    _|_|        _|      _|      \\n"
printf " _|                                                                                \\n"
printf " _|                                                                                \\n"
printf "                                                                                   \\n"
printf " Author: Nobuto NAKAMICHI, Younghwa CHO, Nobuyuki OSHIMA                           \\n"
printf " Date: 25.02.2024                                                                  \\n"
printf " Description:                                                                      \\n"
printf "   This script builds sorce code and creates a new work directory.                 \\n"
printf "\e[0m                                                                              \\n"

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