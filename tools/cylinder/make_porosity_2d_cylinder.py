import matplotlib.pyplot as plt
import numpy as np
import cv2
import seaborn as sns
from scipy.ndimage import convolve
import os

# settings
m,n = 1024, 512
R = 16
thickness = 1.5

def main():

    output_folder = 'data'
    sdf = np.zeros((m,n))

    epsilon_sdf = np.zeros((m,n))
    epsilon_bw_blur = np.zeros((m,n))
    epsilon_gray_blur = np.zeros((m,n))
    epsilon_func = np.zeros((m,n))

    cen_x = (n-1)/2
    cen_y = (n-1)/2

    print(f'center : ({cen_x},{cen_y})')

    for i in range(m):
        for j in range(n):
            dist = abs(np.sqrt((i-cen_x)**2+(j-cen_y)**2))
            sdf[i,j] = dist - R

    bw = sdf_to_bw(sdf)
    gray = sdf_to_gray(sdf)

    os.makedirs(output_folder, exist_ok=True)

    epsilon_func = porosity_func(sdf, thickness=thickness)
    kernel = create_tanh_kernel(thickness=thickness)

    epsilon_gray_blur = convolve(gray, kernel)
    epsilon_bw_blur = convolve(bw, kernel)

    # MSE
    # mse_gray = np.mean(np.abs(epsilon_gray_blur - epsilon_func)**2)
    # mse_bw = np.mean(np.abs(epsilon_bw_blur - epsilon_func)**2)
    # print(f'thickness & MSE (Gray Scale) & MES (Binary)\\\\')
    # print(f'{thickness:.2f} & {mse_gray:.4e} & {mse_bw:.4e}\\\\')

    # epsilon_gray_gauss = cv2.GaussianBlur(gray(sdf), ksize=(ksize, ksize), sigmaX=sigma)
    # epsilon_bw_gauss = cv2.GaussianBlur(gray(sdf), ksize=(ksize, ksize), sigmaX=sigma)


    write_porosity(epsilon_func, f'{output_folder}/func_{thickness}.csv')
    write_porosity(epsilon_gray_blur, f'{output_folder}/gray_blur_{thickness}.csv')
    write_porosity(epsilon_bw_blur, f'{output_folder}/bw_blur_{thickness}.csv')
    # write_porosity(epsilon_gray_gauss, f'{output_folder}/gray_gauss_{thickness}.csv')
    # write_porosity(epsilon_bwgauss, f'{output_folder}/gray_gauss_{thickness}.csv')

def heatmap(data, level, filename='sample'):

    m, n = data.shape
    sz = 16
    print(f'(m, n) = {m,n}')

    sns.heatmap(
        data,
        square=True,
        annot_kws={"size":14},
        # cmap='coolwarm',
        xticklabels=20,
        yticklabels=20,
        )

    # contour lines
    x = np.arange(0.5,m+0.5,1)
    y = np.arange(0.5,n+0.5,1)
    # contour = plt.contour(
    #     x,
    #     y,
    #     data,
    #     levels=level,
    #     colors='black',
    #     linewidths=1
    #     )
    ## contour labelsの追加
    # plt.clabel(contour, inline=True, fontsize=14)

    plt.xlabel('X')
    plt.ylabel('Y')
    # plt.tick_params(axis='both', labelsize=sz)

    plt.savefig(f'./map_{filename}.png')
    plt.clf()

def heatmap_wire(data, filename='sample'):

    m, n = data.shape
    print(f'(m, n) = {m,n}')
    plt.figure(figsize=(8,8))
    sns.heatmap(
        data,
        square=True,
        annot=True,
        annot_kws={"size":14},
        fmt=".3f",
        cbar=False,
        xticklabels=False,
        yticklabels=False,
        cmap="Greys",
        linewidths=1.,
        linecolor="Gray"
        )
    
    plt.savefig(f'./map_{filename}.png')
    plt.clf()


def porosity_func(data, thickness):
    porosity = np.zeros((data.shape))
    porosity = 0.5*np.tanh(data/thickness) + 0.5
    return porosity

def sdf_to_gray(data):
    m, n = data.shape
    data_gray = np.copy(data)
    delta = 1 # px
    for i in range(m):
        for j in range(n):
            val = data[i, j]
            if abs(val) < delta:
                data_gray[i, j] = val/2 + 0.5
            elif val < 0:
                data_gray[i, j] = 0.0
            elif val > 0:
                data_gray[i, j] = 1.0
    return data_gray


def create_tanh_kernel(thickness=2.5):
    size=int(thickness*14) + (1-int(thickness*14)%2)
    x = np.linspace(-size/2, size/2, size)
    y = np.linspace(-size/2, size/2, size)
    x, y = np.meshgrid(x, y)
    kernel = 1 - np.tanh((np.sqrt(x**2 + y**2))/ thickness)
    kernel /= kernel.sum()
    return kernel


def heatmap_3d(data):
    import matplotlib.pyplot as plt
    import numpy as np
    from mpl_toolkits.mplot3d import Axes3D

    x, y = np.meshgrid(np.arange(data.shape[1]), np.arange(data.shape[0]))
    x = x.flatten()
    y = y.flatten()
    z = np.zeros_like(x)

    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.bar3d(x, y, z, dx=1, dy=1, dz=data.flatten())

    ax.set_xlabel('X-axis')
    ax.set_ylabel('Y-axis')
    ax.set_zlabel('Kernel Values')
    # ax.set_title('3D Heatmap using Matplotlib')

    plt.show()


def write_porosity(data,filename='porosity.csv'):
    # Write porosity distribution
    m, n = data.shape
    with open(filename, "w") as file:
        # Write header
        file.write(f"{m},{n},1\n")
        # Write data
        for j in range(1, n+1):
            for i in range(1, m+1):
                file.write(f"{int(i)}, {int(j)}, 1, {(data[i-1, j-1]):.10f}\n")


def sdf_to_bw(sdf):
    m, n = sdf.shape
    bw = np.zeros((m,n))
    for i in range(m):
        for j in range(n):
            if sdf[i, j] > 0:
                bw[i, j] = 1.0
            else:
                bw[i, j] = 0.0
    return bw

    
if __name__ == "__main__":
    main()