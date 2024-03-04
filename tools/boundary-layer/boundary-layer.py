import output_porosity as output
import numpy as np

m, n = 256, 256
thickness = 1.5
data = np.ones((m,n))

for i in range(m):
    for j in range(n):
        dist = abs(128 - j) - 16
        data[i, j] = 0.5*np.tanh(dist/thickness)+0.5

output.write_porosity_2d(data,"porosity.csv")