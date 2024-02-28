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
#   This script helps users by checking available executable files and outputting   #
#   log files.                                                                      #
# _________________________________________________________________________________ #

source monitor.sh
EXE_DIR="../bin/"
LOG_DIR="logs"
STDOUT_FNAME="runlog_$(date "+%Y.%m.%d-%H.%M.%S").txt"

# Create logs directory if it doesn't exist
mkdir -p "$LOG_DIR"

# List available executable files in EXE_DIR
echo "Available executable files:"
counter=0
for file in "$EXE_DIR"*
do
    if [ -x "$file" ]; then
        echo "$counter: ${file##*/}"
        counter=$((counter+1))
    fi
done

# Prompt user to choose an executable file
echo -n "Enter the number of the executable file to run: "
read selection

# Validate user input
if ! [ "$selection" -ge 0 ] || ! [ "$selection" -lt "$counter" ]; then
    echo "\e[1;31mError:\e[0m Invalid selection. Exiting."
    exit 1
fi

# Get the selected executable file
counter=0
selected_file=""
for file in "$EXE_DIR"*
do
    if [ -x "$file" ]; then
        if [ "$counter" -eq "$selection" ]; then
            selected_file="$file"
            break
        fi
        counter=$((counter+1))
    fi
done

# Execute the selected file
if echo "${selected_file##*/}" | grep -q "omp"; then
    echo -n "Enter the number of threads to use for execution:"
    read OMP_NUM_THREADS
    export OMP_NUM_THREADS
    echo " # Number of threads used = $OMP_NUM_THREADS" > "$LOG_DIR/$STDOUT_FNAME"
fi
echo "Running ${selected_file##*/}..."

# Execute the selected file and redirect output to log file
"$selected_file" >> "$LOG_DIR/$STDOUT_FNAME" 2>&1 &

pid=$!

echo "---" >> "$LOG_DIR/process.txt"
echo "The executable file to run: $selected_file" >> "$LOG_DIR/process.txt"
if echo "${selected_file##*/}" | grep -q "omp"; then
    echo "OMP_NUM_THREADS: $OMP_NUM_THREADS" >> "$LOG_DIR/process.txt"
fi
echo ${STDOUT_FNAME} >> "$LOG_DIR/process.txt"
echo "Process PID: $pid" >> "$LOG_DIR/process.txt"