import math
import vtk
import csv
import numpy as np

def main():

    dist_3d_array = process_stl_file_three_axis(
        file_path="sphere.stl",
        bounds_factor=[2.0, 4.0, 1.5, 1.5, 1.5, 1.5],
        grid=60,
        axis=0,
        thickness=1.5
        )
    
    # Save 3D array to CSV file
    save_3d_array_to_csv(
        csv_file_path="output.csv",
        cell_dims=dist_3d_array.shape,
        dist_3d_array=dist_3d_array
        )
    

def ratio_margin_to_bounds(bounds, factor):
    return [bound * factor for bound in bounds]


def ratio_margin_to_bounds_for_three_axis(bounds,bounds_factor):
    
    new_bounds = np.zeros_like(bounds)

    for i in range(6):
        new_bounds[i] = bounds[i] * bounds_factor[i] * 2
    
    return new_bounds


def read_stl_file(file_path):
    reader = vtk.vtkSTLReader()
    reader.SetFileName(file_path)
    reader.Update()
    return reader.GetOutput()


def calculate_pitch_and_mins(bounds, grid, axis):
    pitch = float((bounds[2 * axis + 1] - bounds[2 * axis]) / grid)
    mesh_pitch = [pitch] * 3
    mins = [bound - pitch / 2 for bound in bounds[::2]]
    return pitch, mesh_pitch, mins


def create_mesh_grid_points(cell_dims, mesh_pitch, mins):
    points = vtk.vtkPoints()
    for iz in range(cell_dims[2] + 1):
        for iy in range(cell_dims[1] + 1):
            for ix in range(cell_dims[0] + 1):
                x = ix * mesh_pitch[0] + mins[0]
                y = iy * mesh_pitch[1] + mins[1]
                z = iz * mesh_pitch[2] + mins[2]
                points.InsertNextPoint(x, y, z)
    return points


def convert_to_unstructured_grid(structured_base_mesh):
    append = vtk.vtkAppendFilter()
    append.AddInputData(structured_base_mesh)
    append.Update()
    return append.GetOutput()


def calculate_sdf(poly_data, center_points):
    cell_list = vtk.vtkIdList()
    sdf = vtk.vtkImplicitPolyDataDistance()
    sdf.SetInput(poly_data)
    distance_sdf = vtk.vtkDoubleArray()
    distance_sdf.SetName("sdf")
    num_points = center_points.GetNumberOfPoints()
    point_array = np.array([center_points.GetPoint(idx) for idx in range(num_points)])
    distances = np.array([sdf.FunctionValue(point) for point in point_array])
    for idx, distance in enumerate(distances):
        distance_sdf.InsertNextValue(distance)
        if distance <= 0:
            cell_list.InsertNextId(idx)
    return distance_sdf


def save_3d_array_to_csv(csv_file_path, cell_dims, dist_3d_array):
    with open(csv_file_path, mode='w', newline='') as file:
        writer = csv.writer(file)
        # Write header
        writer.writerow(cell_dims)
        # Write data
        for iz in range(dist_3d_array.shape[2]):
            for iy in range(dist_3d_array.shape[1]):
                for ix in range(dist_3d_array.shape[0]):
                    writer.writerow([ix + 1, iy + 1, iz + 1,\
                                     format(dist_3d_array[ix, iy, iz], '.6E')])


def process_stl_file(file_path, factor, grid, axis, thickness):

    # Read the STL file
    poly_data = read_stl_file(file_path)

    # Add margin to the original size of the STL file
    original_bounds = poly_data.GetBounds()
    expanded_bounds = ratio_margin_to_bounds(original_bounds, factor)
    bounds = expanded_bounds

    # Calculate pitch, mesh_pitch, and mins
    pitch, mesh_pitch, mins = calculate_pitch_and_mins(bounds, grid, axis)

    cell_dims = [math.ceil((bounds[i * 2 + 1] - bounds[i * 2]) / mesh_pitch[i // 2]) for i in range(3)]

    # Creating Mesh Grid Points
    points = create_mesh_grid_points(cell_dims, mesh_pitch, mins)

    structured_base_mesh = vtk.vtkStructuredGrid()
    structured_base_mesh.SetExtent(0, cell_dims[0], 0, cell_dims[1], 0, cell_dims[2])
    structured_base_mesh.SetPoints(points)

    # Convert structured grid data to unstructured grid data
    base_mesh = convert_to_unstructured_grid(structured_base_mesh)

    # Find the coordinates of the center of each Voxel
    cell_centers = vtk.vtkCellCenters()
    cell_centers.SetInputData(base_mesh)
    cell_centers.Update()

    # Create Voxel mesh from STL
    center_points = cell_centers.GetOutput().GetPoints()
    
    distance_sdf = calculate_sdf(poly_data,center_points)
    # Add SDF values
    base_mesh.GetCellData().SetScalars(distance_sdf)

    # Get the cell data array
    cell_data_array = base_mesh.GetCellData().GetArray("sdf")

    dist_3d_array = np.zeros(cell_dims)

    cell_idx = 0
    for iz in range(cell_dims[2]):
        for iy in range(cell_dims[1]):
            for ix in range(cell_dims[0]):
                X = cell_data_array.GetValue(cell_idx) / (thickness * pitch)
                dist_3d_array[ix, iy, iz] = 0.5 * math.tanh(X) + 0.5
                cell_idx += 1

    return dist_3d_array


def process_stl_file_three_axis(file_path, bounds_factor, grid, axis, thickness):

    # Read the STL file
    poly_data = read_stl_file(file_path)

    # Add margin to the original size of the STL file
    original_bounds = poly_data.GetBounds()
    expanded_bounds = ratio_margin_to_bounds_for_three_axis(
        original_bounds,
        bounds_factor
        )
    bounds = expanded_bounds

    # Calculate pitch, mesh_pitch, and mins
    pitch, mesh_pitch, mins = calculate_pitch_and_mins(bounds, grid, axis)

    cell_dims = [math.ceil((bounds[i * 2 + 1] - bounds[i * 2]) / mesh_pitch[i // 2]) for i in range(3)]

    # Creating Mesh Grid Points
    points = create_mesh_grid_points(cell_dims, mesh_pitch, mins)

    structured_base_mesh = vtk.vtkStructuredGrid()
    structured_base_mesh.SetExtent(0, cell_dims[0], 0, cell_dims[1], 0, cell_dims[2])
    structured_base_mesh.SetPoints(points)

    # Convert structured grid data to unstructured grid data
    base_mesh = convert_to_unstructured_grid(structured_base_mesh)

    # Find the coordinates of the center of each Voxel
    cell_centers = vtk.vtkCellCenters()
    cell_centers.SetInputData(base_mesh)
    cell_centers.Update()

    # Create Voxel mesh from STL
    center_points = cell_centers.GetOutput().GetPoints()
    
    distance_sdf = calculate_sdf(poly_data,center_points)
    # Add SDF values
    base_mesh.GetCellData().SetScalars(distance_sdf)

    # Get the cell data array
    cell_data_array = base_mesh.GetCellData().GetArray("sdf")

    dist_3d_array = np.zeros(cell_dims)

    cell_idx = 0
    for iz in range(cell_dims[2]):
        for iy in range(cell_dims[1]):
            for ix in range(cell_dims[0]):
                X = cell_data_array.GetValue(cell_idx) / (thickness * pitch)
                dist_3d_array[ix, iy, iz] = 0.5 * math.tanh(X) + 0.5
                cell_idx += 1

    return dist_3d_array


if __name__ == "__main__":
    main()