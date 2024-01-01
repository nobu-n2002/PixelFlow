#!/bin/bash

# Generate controlDict.txt
cat <<EOF > config/controlDict.txt
&physical
xnue            = 0.025000  ! [m2/s]
xlamda          = 0.000000  ! [m2/s]
density         = 1.000000  ! [kg/m3]
width           = 4.000000  ! [m]
height          = 4.000000  ! [m]
depth           = 4.000000  ! [m]
time            = 10.00000  ! [s]
inlet_velocity  = 1.500000  ! [m/s]
outlet_pressure = 0.000000  ! [gauge]
AoA             = 0.000000  ! [degree]
/
&file_control
istep_out       = 2000
/
&grid_control
istep_max       = 10000
/
&porosity_control
thickness       = 2.0
threshold       = 1.0e-3
/
&calculation_method
nonslip         = .true.  ! .ture.:No-slip cond., .false.:Slip cond.
/
&directory_control
output_folder   = "output"
csv_file        = "data/porosity.csv"
/
EOF

echo "controlDict.txt generated in config directory."
