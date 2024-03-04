import output_porosity as output
import numpy as np
import math

import seaborn as sns
import matplotlib.pyplot as plt


H = 16
wall_thick = 16
step_length = H*10
m, n = H*50 + wall_thick, H*5+wall_thick
cen = n/2

thickness = 1.5
data = np.ones((m,n))

for i in range(m):
    for j in range(n):
        dist = float("inf")

        # upper part
        dist =  min(dist, j - wall_thick/2-cen)

        if i - step_length<= 0:
            dist =  min(dist, j - wall_thick/2 - H - cen)
            if j - wall_thick/2 - H  - cen <= 0:
                dist =  max(dist, j - wall_thick/2 - H  - cen)
            else:
                dist =  min(dist, j - wall_thick/2 - H  - cen)
            
        if j - wall_thick/2 - cen <= 0 and i - step_length <= 0:
            dist = max(dist, -math.sqrt((i-(step_length))*(i-(step_length))+(j-(wall_thick/2+cen))*(j-(wall_thick/2+cen))))

        if j - wall_thick/2 - H - cen <= 0 and j - wall_thick/2 - cen >= 0:
            if i - step_length <= 0:
                dist =  max(dist, i - step_length)
            else:
                dist =  min(dist, i - step_length)

        if j - wall_thick/2 - H - cen >= 0:
            dist = min(dist, math.sqrt((i-step_length)*(i-step_length)+(j-(H+wall_thick/2+cen))*(j-(H+wall_thick/2+cen))))

        # lower part
        dist = max(dist, - (j + wall_thick/2 - cen))

        data[i, j] = min(data[i,j], 0.5*np.tanh(dist/thickness)+0.5)        

output.write_porosity_2d(data,f"../data/porosity_{thickness}_{H}.csv")

sns.heatmap(data, square=True)
plt.show()