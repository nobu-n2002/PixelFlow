def write_porosity_2d (data, filename='porosity.csv'):
    # Write porosity distribution
    m, n = data.shape
    with open(filename, "w") as file:
        # Write header
        file.write(f"{m},{n},1\n")
        # Write data
        for j in range(1, n+1):
            for i in range(1, m+1):
                file.write(f"{int(i)}, {int(j)}, 1, {(data[i-1, j-1]):.10f}\n")


def write_porosity_3d (data,filename='porosity.csv'):
    # Write porosity distribution
    m, n, l = data.shape
    with open(filename, "w") as file:
        # Write header
        file.write(f"{m},{n},{l}\n")
        # Write data
        for k in range(1, l+1):
            for j in range(1, n+1):
                for i in range(1, m+1):
                    file.write(f"{int(i)}, {int(j)}, {int(k)}, {(data[i-1, j-1, k-1]):.10f}\n")


