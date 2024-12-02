![GitHub code size in bytes](https://img.shields.io/github/languages/code-size/nobu-n2002/PixelFlow)
![Static Badge](https://img.shields.io/badge/Ubuntu_v20.04-supported-blue)


# PixelFlow

## Overview

This implementation is based on the [papers](#References).

This code deals with incompressible fluids and implements the collocated grid Marker-and-Cell (MAC) method using a regular orthogonal grid. It discretizes space using second-order central differencing and time using a first-order explicit Euler method. For pressure calculation, it uses the Red-Black Successive Over-Relaxation (SOR) method.

## Table of Contents

- [PixelFlow](#pixelflow)
  - [Overview](#overview)
  - [Table of Contents](#table-of-contents)
  - [Requirements](#requirements)
    - [Prerequisites](#prerequisites)
  - [Installation](#installation)
  - [Getting Started](#getting-started)
    - [Setting up the Environment](#setting-up-the-environment)
    - [Editing Configuration](#editing-configuration)
    - [Running Simulations](#running-simulations)
  - [Configuring Simulations](#configuring-simulations)
    - [\&physical Section](#physical-section)
    - [\&file\_control Section](#file_control-section)
    - [\&grid\_control Section](#grid_control-section)
    - [\&porosity\_control Section](#porosity_control-section)
    - [\&calculation\_method Section](#calculation_method-section)
    - [\&directory\_control Section](#directory_control-section)
  - [Running a Test Case](#running-a-test-case)
  - [References](#references)


## Requirements

To build and run this project, you will need the following tools and libraries installed on a Linux environment:

- **CMake**: A cross-platform build system generator.
- **Make**: A build automation tool.
- **GCC**: The GNU Compiler Collection, including the C compiler.
- **G++**: The GNU C++ Compiler.
- **GFortran**: The GNU Fortran Compiler.

### Prerequisites

Make sure you have the required tools and libraries installed. You can install them using the package manager of your Linux distribution. For example, on Debian-based systems like Ubuntu, you can use `apt-get`:

```bash
sudo apt-get update
sudo apt-get install -y cmake make gcc g++ gfortran
```

## Installation

1. Open a terminal and navigate to the desired working directory.

2. Clone the repository from GitHub:

    ```bash
    git clone https://github.com/nobu-n2002/PixelFlow.git
    ```

3. Move into the application directory:

    ```bash
    cd PixelFlow
    ```

4. Build all sources:

    ```bash
    make build
    ```

5. Create working directory:

    ```bash
    make project f=<your_project_name>
    ```

    The usage of `make` is as follows:

    ```makefile
    Usage: make [OPTIONS]
    options
    build                    Build all source codes
    project f=<folder_name>  Create a new folder with the specified name
    help                     Display Usage
    ```

6. Check the `bin` Folder:
   - Upon a successful build, executable files are generated in the `bin` folder.
   - Open a terminal or command prompt and run the following command to verify the presence of the generated executable file(s) in the `bin` folder.

     ```bash
     ls bin
     ```

   - Ensure that the application's executable file is present in the `bin` folder.

7. Grant execution permissions to the executable files in the `bin` folder:

    ```bash
    chmod +x bin/*
    ```

8. Verify the `config` Folder:
   - Proper configuration is crucial for the application. During the build, a configuration file named `controlDict.txt` should be created in the `config` folder.
   - Open a terminal or command prompt and run the following command to confirm the existence of the `controlDict.txt` file in the `config` folder.

     ```bash
     cd projects/<folder_name>
     ls config
     ```

   - Confirm that the `controlDict.txt` file is present in the `config` folder.

9. Follow the steps in the [Getting Started](#getting-started) section to set up the environment and configure the simulation.

10. Proceed to the [Running Simulations](#running-simulations) section to execute the simulation.

## Getting Started

### Setting up the Environment

1. Place the solid boundary information file (`.csv`) in the `data/` directory.

### Editing Configuration

1. Open the `config/controlDict.txt` file in a text editor. For detailed information on each parameter and how to configure your simulations, refer to the [Configuring Simulations](#configuring-simulations) section.

2. Edit the variables according to your simulation requirements.

    ```plaintext
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
    ...
   !********************************************
   &solver_control
   !--- SOR max iteration steps
   iter_max        = 100
   !--- SOR reluxation factor (1<w<2)
   relux_factor    = 1.700000
   /
   !********************************************
    ```

3. Save the changes and close the file.

### Running Simulations

1. Run the simulation script:

    ```bash
    make run
    ```

   You will be prompted to select an executable file. If you are using OpenMP parallelization code, please set parallel threads in `config/omp_config.conf`.

   ```bash
   Available executable files:
   0: ibm2_drag_omp
   1: ibm2_omp
   2: ibm3_air_condition_omp
   3: ibm3_omp
   Using OMP_NUM_THREADS = ${$OMP_NUM_THREADS}
   Running ibm2_omp...
   ```

2. The processes are output to `runlog_*.log` files inside the `logs/` folder.
   Additionally, execution information is output to `process.log` inside the `logs/` folder.
3. The output in the `{output_folder}/` directory.

4. If you need to forcibly terminate the process during computation, please enter the following command. You will receive the process ID and a confirmation message. Input 'yes' followed by pressing enter to forcibly terminate the running process.

    ```bash
    make quit
    ```

## Configuring Simulations

The `config/controlDict.txt` file contains parameters that define the properties of the simulation. Here's a breakdown of the important parameters:

### &physical Section

- **`xnue`**: Kinematic viscosity of the fluid, $\nu = (\mu/\rho$) [m^2/s].
- **`xlambda`**: Second Kinematic viscosity of the fluid, ($\lambda/\rho$) [m^2/s].
- **`density`**: Density of the fluid, $\rho$ [kg/m^3].
- **`width`**, **`height`**, **`depth`**: The simulation domain, $X, Y, Z$ [m].
- **`time`**: Total simulation time, $t$ [s].
- **`inlet_velocity`**: Inlet velocity of the fluid, $U_{in}$ [m/s].
- **`outlet_pressure`**: Outlet pressure, $P_{out}$ [Pa].
- **`AoA`**: Angle of Attack, $\alpha$ [degree].

### &file_control Section

- **`istep_out`**: Output interval for saving results, specified as the number of time steps.

### &grid_control Section

- **`istep_max`**: Maximum number of time steps for the simulation.

### &porosity_control Section

- **`thickness`**: Thickness of boundary region, $\Delta^*(=\Delta/dx$).

### &calculation_method Section

- **`nonslip`**: Specify `.true.` for no-slip conditions or `.false.` for slip conditions.

### &directory_control Section

- **`output_folder`**: Name of the folder where simulation results will be stored.
- **`csv_file`**: Path to the solid boundary information file in CSV format.

Adjust these parameters according to your simulation requirements. The `output_folder` will be created to store the simulation results, and the `csv_file` should point to the CSV file containing solid boundary information.

## Running a Test Case

- [2D circler cylinder case](doc/test-case/cylinder-2d)

## References

[1] Oshima, N., A novel approach for wall-boundary immersed flow simulation (proposal of modified Navier-Stokes equation), Journal of Fluid Science and Technology. Vol.18, No.4 (2023), 
https://doi.org/10.1299/jfst.2023jfst0034

[2] Oshima, N., A novel approach for wall-boundary immersed flow simulation (part 2: modeling of wall shear stress), Journal of Fluid Science and Technology. Vol.19, No.3 (2024), https://doi.org/10.1299/jfst.2024jfst0026

[3] Oshima, N., Program for flow simulation immersing wall boundary, Hokkaido university collection of scholarly and academic papers, http://hdl.handle.net/2115/89344

[4] Nakamichi, N., Cho, Y. and Oshima, N., Image-data-driven simulation of fluid dynamics (proposal and evaluation), Mechanical Engineering Journal, Advance online publication, https://doi.org/10.1299/mej.24-00196
