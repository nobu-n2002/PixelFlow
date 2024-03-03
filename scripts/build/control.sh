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
#   This script extracts the necessary controlDict.txt file and the executable      #
#   script, serving as the execution tool, into a new project folder.               #
# _________________________________________________________________________________ #

# Generate controlDict.txt
cat <<EOF > config/controlDict.txt
!********************************************
&physical
!--- Kinematic viscosity coefficient [m2/s]
xnue = 0.001000
!--- Second viscosity coefficient [m2/s]
xlambda = 0.000000
!--- Fluid density [kg/m3]
density = 1.000000
!--- Domein [m]
width   = 1.0000000
height  = 1.0000000
depth   = 1.0000000
!--- Simulation time [s]
time = 1.000000 
!--- Inlet velosity [m/s]
inlet_velocity = 0.100000
!--- Outlet pressure [Pa]
outlet_pressure = 0.000000
!--- Angle of Attack [degree]
AoA = 0.000000
/
!********************************************
&file_control
istep_out = 1001
/
!********************************************
&grid_control
istep_max = 100
/
!********************************************
&porosity_control
!--- Delta/dx
thickness = 1.500000
!--- munimum porosity value
threshold = 1.0e-6
!--- object radius [m]
radius = 0.100000
!--- relative position from the coordinate center
center_x = 0.500000
center_y = 0.500000
center_z = 0.500000
/
!********************************************
&calculation_method
!--- .ture.: No-slip cond. | .false.: Slip cond.
nonslip = .true.
/
!********************************************
&directory_control
!--- file/folder name (within 50 characters)
output_folder   = "output"
csv_file        = "data/porosity.csv"
/
!********************************************
&solver_control
!--- SOR max iteration steps
iter_max        = 100
!--- SOR reluxation factor (1<w<2)
relux_factor    = 1.700000
/
!********************************************
EOF

echo "controlDict.txt generated in config directory."

cp ../scripts/run.sh run.sh
echo "Copied the execution script."

cp ../scripts/quit.sh quit.sh
echo "Copied the termination script."