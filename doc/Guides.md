# PixelFlow

## Overview

This is an implementation based on the paper "A Novel approach for wall-boundary immersed flow simulation: proposal of modified Navier-Stokes equation" by Nobuyuki OSHIMA, published in the Mechanical Engineering Journal, Volume 18, Number 4 (2023).

## Citing A Novel approach for wall-boundary immersed flow simulation

For comprehensive insights into the proposed methodology and findings presented in our work, please consider referencing the [paper](https://doi.org/10.1299/jfst.2023jfst0034):

```bibtex
@InProceedings{Oshima_2023_MEJ,
  author    = {Nobuyuki OSHIMA},
  title     = {A Novel Approach for Wall-boundary Immersed Flow Simulation: Proposal of Modified Navier-Stokes Equation},
  booktitle = {Mechanical Engineering Journal},
  volume    = {18},
  number    = {4},
  year      = {2023},
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
  - [References](#references)
  - [Example](#example)
    - [Running a Test Case](#running-a-test-case)
    - [Expected Output](#expected-output)

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
    chmod +x init.sh run.sh
    ```

5. Run the initialization script:

    ```bash
    sh init.sh
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
    chmod +x bin/ibm2 bin/ibm3
    ```

8. Verify the `config` Folder:
   - Proper configuration is crucial for the application. During the build, a configuration file named `controlDict.txt` should be created in the `config` folder.
   - Open a terminal or command prompt and run the following command to confirm the existence of the `controlDict.txt` file in the `config` folder.

     ```bash
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
    &physical
    xnue = 0.025000
    # ... other parameters ...
    /

    &file_control
    istep_out = 10001
    /
    # ... other sections ...
    ```

3. Save the changes and close the file.

### Running Simulations

1. Open the `run.sh` script in a text editor.

2. Locate the `DIMENSION` variable within the script and set it to either 2 or 3.

    ```bash
    # init.sh

    # Set the dimension to 2 or 3
    DIMENSION=2
    ```

3. Run the simulation script:

    ```bash
    sh run.sh
    ```

4. View the progress in the `logs/` directory, and the output in the `{output_folder}/` directory.
5. If you need to forcibly stop the computation midway, you can check the running processes by entering ps in the terminal. For example, the output might look like this:

    ```bash
    ps
    ```

   ```bash
   Copy code
     PID TTY          TIME CMD
   16149 pts/6    00:00:00 bash
   18326 pts/6    00:38:34 ibm2
   18931 pts/6    00:00:00 ps
   ```

   To stop the execution of ibm2, input kill {PID} in the terminal. This will terminate the running process of ibm2.

    ```bash
    kill 18326
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

## References

[1] Oshima.N, A Novel approach for wall-boundary immersed flow simulation: proposal of modified Navier-Stokes equation, Mechanical Engineering Journal. Vol.18, No.4 (2023)

[2] 大島, 流れの数値解析:固体境界が埋め込まれた改良ナビエ・ストークス方程式の解法, 北海道大学学術成果コレクション(HUBCAP), 資源タイプsoftware (2023), [Link](http://hdl.handle.net/2115/89344)

## Example

### Running a Test Case

To ensure that the application is set up correctly, you can run a provided test case located in the `test/` folder. Follow these steps:

1. Open a terminal and navigate to the application directory:

    ```bash
    cd PixelFlow
    ```

2. Run the provided test script and grant execution permissions to `init.sh` and `run.sh`:

    ```bash
    cd test
    ```

    ```bash
    chmod +x init.sh run.sh
    ```

3. Run the initialization script:

    ```bash
    sh init.sh
    ```

4. Confirm the creation of the executable files in the `bin` folder and the `controlDict.txt` file in the `config` folder:

    ```bash
    ls bin/ibm2 config/controlDict.txt
    ```

5. Grant execution permissions to the executable files in the `bin` folder:

    ```bash
    chmod +x bin/ibm2
    ```

6. Run the simulation script:

    ```bash
    sh run.sh
    ```

7. Monitor the progress in the `logs/` directory and check for successful execution in the `{output_folder}/` directory.
8. Upon successful completion of the 2D test case, your directory structure should resemble the following:

   ```plaintext
   .
   ├── bin
   │   └── ibm2
   ├── config
   │   └── controlDict.txt
   ├── data
   │   ├── circle.csv
   │   └── stanford_dragon.csv
   ├── data.zip
   ├── etc
   │   ├── divergent.dat
   │   ├── grid.dat
   │   ├── solution_uvp.dat
   │   └── surface_profile.dat
   ├── init.sh
   ├── logs
   │   └── runlog_*.txt
   ├── run.sh
   ├── scripts
   │   ├── SQUID
   │   │   ├── qcompile.sh
   │   │   └── qrun.sh
   │   ├── build2d.sh
   │   ├── build3d.sh
   │   ├── cleanAll.sh
   │   ├── test2d.sh
   │   └── test3d.sh
   ├── src
   │   ├── ibm_2d.f90
   │   └── ibm_3d.f90
   └── test2d
       └── output_paraview.vtk
   ```

9. To confirm the successful execution of the 2D test case, follow these steps:

   1. Locate the `output_paraview.vtk` file in the `test2d` directory.

   2. Open the file using visualization software such as Paraview.

      - If you can visualize the results in Paraview without any issues, it indicates the completion of the test case execution.

   Congratulations! You have successfully completed the 2D test case. If you encounter any difficulties or have further questions, please refer to the provided documentation or seek assistance from our support channels.

   Thank you for your efforts!

10. To perform computations for a 3D test case, set DIMENSION=3 in both init.sh and run.sh, and execute the following commands in sequence:

    ```bash
    # test/init.sh

    # Set the DIMENSION variable to 2 for a 2D cylindrical case
    # Set the DIMENSION variable to 3 for a 3D Stanford Dragon case
    DIMENSION=3
    ```

    ```bash
    # test/run.sh

    # Set the DIMENSION variable to 2 for a 2D cylindrical case
    # Set the DIMENSION variable to 3 for a 3D Stanford Dragon case
    DIMENSION=3
    ```

    ```bash
    sh init.sh
    ```

    ```bash
    sh run.sh
    ```

### Expected Output

Upon successful execution of the test case, you should observe the following:

- A log file in the `logs/` directory (e.g., `logs/runlog_test.txt`).
- Output files in the `test2d` or `test3d` directory.

Feel free to explore the contents of the log file and output folder to verify that the simulation ran as expected.
