# PixelFlow

## Overview

This is an implementation based on the paper "A Novel approach for wall-boundary immersed flow simulation: proposal of modified Navier-Stokes equation" by Nobuyuki OSHIMA, published in the Mechanical Engineering Journal, Volume 18, Number 4 (2023).

This code deals with incompressible fluids and implements the collocated grid MAC method using a regular orthogonal grid. It discretizes space using second-order central differencing and time using a first-order explicit Euler method. For pressure calculation, it implements the Red-Black SOR method.

## Citing A Novel approach for wall-boundary immersed flow simulation

For comprehensive insights into the proposed methodology and findings presented in our work, please consider referencing the [paper](https://doi.org/10.1299/jfst.2023jfst0034):

```bibtex
@article{Oshima_2023_MEJ,
  title={A novel approach for wall-boundary immersed flow simulation},
  author={Nobuyuki OSHIMA},
  journal={Journal of Fluid Science and Technology},
  volume={18},
  number={4},
  pages={23-00192},
  year={2023},
  doi={10.1299/jfst.2023jfst0034}
}
```

## Table of Contents

- [PixelFlow](#pixelflow)
  - [Overview](#overview)
  - [Citing A Novel approach for wall-boundary immersed flow simulation](#citing-a-novel-approach-for-wall-boundary-immersed-flow-simulation)
  - [Table of Contents](#table-of-contents)
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

4. Make the initialization script executable:

    ```bash
    chmod +x project.sh
    ```

5. Run the initialization script:

    ```bash
    sh project.sh -b -f your_project
    ```

    The usage of project.sh is as follows:

    ```bash
   Usage: project.sh [OPTIONS] <folder_name>
   Options:
     -b                 Build all source codes
     -f <folder_name>   Create a new folder with the specified name
     -h                 Display Usage
    ```

    If you only want to build the source code, please execute it with only the `-b` option as follows.

    ```bash
    sh project.sh -b
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
     cd your_project
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
    sh run.sh
    ```

   You will be prompted to select an executable file. If you are using OpenMP parallelization code, please enter the number of OpenMP parallel threads.

   ```bash
   Available executable files:
   0: ibm2_drag_omp
   1: ibm2_omp
   2: ibm3_air_condition_omp
   3: ibm3_omp
   Enter the number of the executable file to run: 1
   Enter the number of threads to use for execution:3
   Running ibm2_omp...
   ```

2. The processes are output to `runlog_*.txt` files inside the `logs/` folder.
   Additionally, execution information is output to `process.txt` inside the `logs/` folder.
3. The output in the `{output_folder}/` directory.
4. If you need to forcibly stop the computation midway, you can check the running processes by entering ps in the terminal. For example, the output might look like this:

    ```bash
    sh quit.sh
    ```

## Configuring Simulations

The `config/controlDict.txt` file contains parameters that define the properties of the simulation. Here's a breakdown of the important parameters:

### &physical Section

- **`xnue`**: Kinematic viscosity of the fluid in [m^2/s].
- **`xlamda`**: Second Kinematic viscosity of the fluid in[m^2/s].
- **`density`**: Density of the fluid in [kg/m^3].
- **`width`**, **`height`**, **`depth`**: The simulation domain in [m].
- **`time`**: Total simulation time in [s].
- **`inlet_velocity`**: Inlet velocity of the fluid in [m/s].
- **`outlet_pressure`**: Outlet pressure in [gauge].
- **`AoA`**: Angle of Attack in [degree].

### &file_control Section

- **`istep_out`**: Output interval for saving results, specified as the number of time steps.

### &grid_control Section

- **`istep_max`**: Maximum number of time steps for the simulation.

### &porosity_control Section

- **`thickness`**: Thickness of boundary region ($\Delta/dx$).

### &calculation_method Section

- **`nonslip`**: Specify `.true.` for no-slip conditions or `.false.` for slip conditions.

### &directory_control Section

- **`output_folder`**: Name of the folder where simulation results will be stored.
- **`csv_file`**: Path to the solid boundary information file in CSV format.

Adjust these parameters according to your simulation requirements. The `output_folder` will be created to store the simulation results, and the `csv_file` should point to the CSV file containing solid boundary information.

## Running a Test Case

- [2D circler cylinder case](doc/test-case/cylinder-2d)

## References

[1] Oshima.N, A Novel approach for wall-boundary immersed flow simulation: proposal of modified Navier-Stokes equation, Mechanical Engineering Journal. Vol.18, No.4 (2023)

[2] 大島, 流れの数値解析:固体境界が埋め込まれた改良ナビエ・ストークス方程式の解法, 北海道大学学術成果コレクション(HUBCAP), 資源タイプsoftware (2023), [Link](http://hdl.handle.net/2115/89344)
