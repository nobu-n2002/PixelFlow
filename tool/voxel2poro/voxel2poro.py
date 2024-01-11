from PIL import Image
import numpy as np
import os
import pyvista as pv
from tqdm import tqdm
from scipy.ndimage import convolve
from pyevtk.hl import imageToVTK


input_folder = 'sample'
output_folder = f'output_view_{input_folder}'
dim = 64
buff = 23
thickness = 1.5
output_vtk = True
output_porosity = True
output_file_path = f'{output_folder}/porosity_{input_folder}.csv'


def load_bitmap_image(file_path):
    # Read Bitmap image as PIL Image Object
    image = Image.open(file_path)
    # Convert to NumPy
    image_array = np.array(image)
    # Change values into binary
    image_array[image_array == 0.] = 1.0
    image_array[image_array == 128.] = 0.0
    image_array[image_array == 255.] = 0.0
    return image_array

def view_threshold(array, output_folder, filename, threshold):
    # PyVista：StructuredGrid
    grid = pv.ImageData()
    grid.dimensions = array.shape
    grid.origin = (0, 0, 0)
    grid.spacing = (1, 1, 1)
    # Data
    grid.point_data['Values'] = array.flatten(order='F').astype(np.float32)
    # Plotter
    plotter = pv.Plotter(off_screen=True)
    plotter.add_mesh(grid.threshold(threshold), show_edges=False)
    plotter.camera_position = [3, 3.5, 1]
    plotter.remove_scalar_bar()
    plotter.add_axes()
    # Create output folder
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    # Save the figure as HTML file
    output_path_html = os.path.join(output_folder, f'view_{filename}.html')
    plotter.export_html(output_path_html)
    # Save the figure as PNG file
    output_path_png = os.path.join(output_folder, f'view_{filename}.png')
    plotter.screenshot(output_path_png)


def view_contour(array, output_folder, filename, isovalue):
    # PyVista：StructuredGrid
    grid = pv.ImageData()
    grid.dimensions = array.shape
    grid.origin = (0, 0, 0)
    grid.spacing = (1, 1, 1)
    # Data
    grid.point_data['Values'] = array.flatten(order='F').astype(np.float32)
    # Plotter
    plotter = pv.Plotter(off_screen=True)
    plotter.add_mesh(grid.contour(isosurfaces=[isovalue]), show_edges=False)
    plotter.camera_position = [3, 3.5, 1]
    plotter.remove_scalar_bar()
    plotter.add_axes()

    # Create output folder
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    # Save the figure as HTML file
    output_path_html = os.path.join(output_folder, f'view_{filename}.html')
    plotter.export_html(output_path_html)
    # Save the figure as PNG file
    output_path_png = os.path.join(output_folder, f'view_{filename}.png')
    plotter.screenshot(output_path_png)


def view_threshold_all(array, output_folder, filename, threshold):
    # PyVista：StructuredGrid
    grid = pv.ImageData()
    grid.dimensions = array.shape
    grid.origin = (0, 0, 0)
    grid.spacing = (1, 1, 1)
    grid.point_data['Values'] = array.flatten(order='F').astype(np.float32)

    # Data
    mesh = grid.threshold(threshold)

    # Subplot Settings
    plot_shape = (2, 3)
    plotter = pv.Plotter(shape=plot_shape, off_screen=True)
    plane = ['xy', 'xz', 'zy', 'yx', 'zx', 'yz']
    roll = [0, 0, 90, 0, 180, -90]

    for i in range(int(plot_shape[0]*plot_shape[1])):
        column = int(i % plot_shape[1])
        row = 0 if i < plot_shape[1] else 1
        plotter.subplot(row, column)
        plotter.add_mesh(mesh, show_edges=False)
        plotter.camera_position = plane[i]
        plotter.camera.roll = roll[i]
        plotter.remove_scalar_bar()
        plotter.add_axes()      

    # Create output folder
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    # Save the figure as PNG file
    output_path = os.path.join(output_folder, f'view_threshold_all_{filename}.png')
    plotter.screenshot(output_path)


def view_contour_all(array, output_folder, filename, isovalue):
    # PyVistaのStructuredGridに変換
    grid = pv.ImageData()
    grid.dimensions = array.shape
    grid.origin = (0, 0, 0)
    grid.spacing = (1, 1, 1)
    grid.point_data['Values'] = array.flatten(order='F').astype(np.float32)

    mesh = grid.contour(isosurfaces=[isovalue])

    # Subplot Settings
    plot_shape = (2, 3)
    plotter = pv.Plotter(shape=plot_shape, off_screen=True)
    plane = ['xy', 'xz', 'zy', 'yx', 'zx', 'yz']
    roll = [0, 0, 90, 0, 180, -90]

    for i in range(int(plot_shape[0]*plot_shape[1])):
        column = int(i % plot_shape[1])
        row = 0 if i < plot_shape[1] else 1
        plotter.subplot(row, column)
        plotter.add_mesh(mesh, show_edges=False)
        plotter.camera_position = plane[i]
        plotter.camera.roll = roll[i]
        plotter.remove_scalar_bar()
        plotter.add_axes() 

    # Create output folder
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    # Save the figure as PNG file
    output_path = os.path.join(output_folder, f'view_contour_all_{filename}.png')
    plotter.screenshot(output_path)


def create_tanh_kernel(thickness=2.0):
    size=int(thickness*4) + (1-int(thickness*4)%2)
    x = np.linspace(-size/2, size/2, size)
    y = np.linspace(-size/2, size/2, size)
    z = np.linspace(-size/2, size/2, size)
    x, y, z = np.meshgrid(x, y, z)
    kernel = 1 - np.tanh((np.sqrt(x**2 + y**2 + z**2))/ thickness)
    kernel /= kernel.sum()
    return kernel


def write_porosity(data,filename='porosity.csv'):
    # Write porosity distribution
    m, n, l = data.shape
    with open(filename, "w") as file:
        # Write header
        file.write(f"{m},{n},{l}\n")
        # Write data
        for k in tqdm(range(1, l+1)):
            for j in range(1, n+1):
                for i in range(1, m+1):
                    file.write(f"{int(i)}, {int(j)}, {int(k)}, {(data[i-1, j-1, k-1]):.10f}\n")


def main():
    array_3d = np.ones((dim+2*space, dim+2*space, dim+2*space), dtype=np.float32)

    for i in range(dim):
        # Load bitmap image by specifying its path
        bitmap_path = f'{input_folder}/img{i:05d}.bmp'
        bitmap_array = load_bitmap_image(bitmap_path)
        array_3d[space:space+dim, space:space+dim, space+i] = bitmap_array[:, :]

    # Visualization    
    view_threshold(array_3d, output_folder, "binary", [-0.1,0.1])
    view_threshold_all(array_3d, output_folder, "binary", [-0.1,0.1])

    print("Calculating porosity value using Tanh Filter...")
    kernel = create_tanh_kernel(thickness=thickness)
    porosity = convolve(array_3d, kernel, mode='nearest', cval=1.000000)
    print(f'porosity shape: {porosity.shape}')
    print("Processing complete.")

    if output_vtk:
        print("Saving Voxel Data to VTK file...")
        imageToVTK(f'{output_folder}/point_data_{input_folder}',pointData={'porosity':porosity})
        print("Processing complete.")

    if output_porosity:
        print('Saving porosity Data to CSV file')
        write_porosity(porosity,filename=output_file_path)
        print(f"porosity data has been saved to {output_file_path}.")

    # Visualization
    view_contour(porosity, output_folder, "porosity0119", 0.119)
    view_contour(porosity, output_folder, "porosity0500", 0.500)
    view_contour(porosity, output_folder, "porosity0881", 0.881)
    view_contour_all(porosity, output_folder, "porosity0119", 0.119)
    view_contour_all(porosity, output_folder, "porosity0500", 0.500)
    view_contour_all(porosity, output_folder, "porosity0881", 0.881)

if __name__ == "__main__":
    main()
