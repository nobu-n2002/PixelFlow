import numpy as np
from PIL import Image
import os
import glob

def create_voxel_data(size, radius):
    # Sphere case
    data = np.zeros((size, size, size), dtype=np.uint8)
    center = size // 2 - 0.5

    for x in range(size):
        for y in range(size):
            for z in range(size):
                distance = np.sqrt((x - center)**2 + (y - center)**2 + (z - center)**2)
                if distance <= radius:
                    data[x, y, z] = 255.0  # 白色

    return data

def save_bitmap(data, output_path):
    img = Image.fromarray(data, 'L')
    img.save(output_path)

def main():
    voxel_size = 64
    sphere_radius = 16  # 半径
    num_spheres = 64

    # Create output folder
    output_folder = 'sample'
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    files_to_delete = glob.glob(f"{output_folder}/*.bmp")
    for file_path in files_to_delete:
        os.remove(file_path)

    for i in range(num_spheres):
        voxel_data = create_voxel_data(voxel_size, sphere_radius)
        output_path = f"{output_folder}/img{i:05d}.bmp"
        save_bitmap(voxel_data[i, :, :], output_path)  # 32x32のスライスを使用
        print(f"Saved {output_path}")

if __name__ == "__main__":
    main()
