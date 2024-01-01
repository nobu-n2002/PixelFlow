#!/bin/bash

# Generate controlDict.txt
cat <<EOF > config/controlDict.txt
&physical
xnue            = 0.025000  ! [m2/s]
xlamda          = 0.000000  ! [m2/s]
density         = 1.000000  ! [kg/m3]
width           = 5.120000  ! [m]
height          = 3.620000  ! [m]
depth           = 2.300000  ! [m]
time            = 6.250000  ! [s]
inlet_velocity  = 1.500000  ! [m/s]
outlet_pressure = 0.000000  ! [gauge]
AoA             = 0.000000  ! [degree]
/
&file_control
istep_out       = 1000
/
&grid_control
istep_max       = 5000
/
&porosity_control
thickness       = 2.0
/
&calculation_method
nonslip         = .true.  ! .ture.:No-slip cond., .false.:Slip cond.
/
&directory_control
output_folder   = "test3d"
csv_file        = "data/stanford_dragon.csv"
/
EOF

echo "controlDict.txt generated in config directory."
