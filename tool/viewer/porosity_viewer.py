from pyevtk.hl import imageToVTK
import pandas as pd
import numpy as np
import os


def main():
    porosity_array = read_porosity(input_porosity_file='../stl2poro/output.csv')
    porosity2vtk(
        input_porosity_array=porosity_array,
        output_folder='output',
        output_file_name='porosity',
    )


def read_porosity(input_porosity_file):
    # Read porosity data from CSV file
    grid = pd.read_csv(input_porosity_file, header=None, nrows=1)

    m, n, l = grid.iloc[0, :3].astype(int)
    print(f'(m, n, l) = {m,n,l}')

    porosity_data = pd.read_csv(input_porosity_file, skiprows=1, usecols=[0, 1, 2, 3], header=None)

    # Obtain the maximum values for x, y

    porosity = np.ones((m, n, l))


    # Extract porosity values
    for i in range(int(m*n*l)):
        x = porosity_data[0][i]
        y = porosity_data[1][i]
        z = porosity_data[2][i]
        val = porosity_data[3][i]
        porosity[x-1, y-1, z-1] = val

    return porosity


def porosity2vtk(input_porosity_array, output_folder, output_file_name):
    os.makedirs(output_folder, exist_ok=True)
    print("Saving Voxel Data to VTK file...")
    imageToVTK(
        f'{output_folder}/point_data_{output_file_name}',
        pointData={'porosity':input_porosity_array}
        )
    print("Processing complete.")


if __name__ == '__main__':
    main()