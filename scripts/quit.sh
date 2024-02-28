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

for i in $pid
do
    # check processes
    if ps $i > /dev/null; then
        echo "Process with PID $i is running."
        # Prompt for user input with a warning message
        echo -n "Are you sure you want to terminate? (yes/no): "
        # Read user input
        read input
        # If user input is "yes"
        if [ "$input" = "yes" ] || [ "$input" = "y" ]; then
            # Terminate the process using the extracted PID
            kill $i
            echo "Process with PID $i has been terminated."
            exit 0
        # If user input is "no"
        elif [ "$input" = "no" ] || [ "$input" = "n" ]; then
            echo "Not terminating."
            exit 1
        # If user input is neither "yes" nor "no"
        else
            echo "Invalid input. Please enter 'yes' or 'no'."
            exit 1
        fi
    else
        echo "Process with PID $i is not running."
    fi
done