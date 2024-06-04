import output_porosity as output
import numpy as np
import math
import seaborn as sns
import matplotlib.pyplot as plt

def main():
    
    H = 16
    h = 128
    m, n = H*8 + h, H*4+h
    cen = n/2

    thickness = 1.5
    data = np.ones((m,n))

    for i in range(m):
        for j in range(n):
            dist = float("inf")

            # upper part
            dist =  min(dist, - (j - h/2-cen))

            # lower part
            dist = min(dist, (j + h/2 - cen))

            data[i, j] = (np.tanh(dist/thickness))**2        

    output.write_porosity_2d(data,f"porosity_{thickness}_{H}.csv")

    sns.heatmap(data, square=True)
    plt.show()
    
if __name__ == "__main__":
    main()