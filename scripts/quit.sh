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
#   This script checks processes and terminates them.                               #
# _________________________________________________________________________________ #

# Extract PID from "process.txt" file
pid=$(grep -oP 'Process PID: \K\d+' logs/process.txt)

echo "PROCESS PID :STATUS"
echo "---------------------"
for i in $pid
do
    # check processes
    if ps $i > /dev/null; then
        echo "PID $i :RUN"
        while true; do
            # Prompt for user input with a warning message
            echo -n "--- Are you sure you want to terminate $i? (yes/no): "
            # Read user input
            read input
            case "$input" in
                [yY]|[yY][eE][sS])
                    # Terminate the process using the extracted PID
                    kill $i
                    echo "Process with PID $i has been terminated."
                    break
                    ;;
                [nN]|[nN][oO])
                    echo "Terminating canceled."
                    break
                    ;;
                *)
                    echo "Invalid input. Please enter 'yes' or 'no'."
                    ;;
            esac
        done
    else
        echo "PID $i :END"
    fi
done