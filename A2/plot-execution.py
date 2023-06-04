import pandas as pd
import numpy as np
import glob
import matplotlib.pyplot as plt


def get_arguments_from_filename(filename: str) -> str:
    array = filename.replace(".csv", "").split("-")
    n = array[3]
    m = array[4]
    max_iterations = array[5]
    epsilon = array[6]

    return f"n={n}, m={m}, max_iter={max_iterations}, \u03B5={epsilon}"


omp_execution_time_files = glob.glob("execution-time-omp-*")
mpi_execution_time_files = glob.glob("execution-time-mpi-*")

plt.figure(figsize=(15, 8))

plt.subplot(2, 1, 1)
plt.grid()
plt.title("OMP speedup")
plt.xlabel("number of threads")
plt.ylabel("speedup")
plt.xticks(np.arange(1, 33))
plt.yticks(np.arange(1, 20))

for filename in omp_execution_time_files:
    df = pd.read_csv(filename)
    plt.plot(df["num_thread"], df["execution-time"].iloc[0] / df["execution-time"], label=get_arguments_from_filename(filename))

plt.legend()
plt.subplot(2, 1, 2)
plt.title("MPI speedup")
plt.xlabel("configuration (nodes, threads)")
plt.ylabel("speedup")
plt.grid()
for filename in mpi_execution_time_files:
    df = pd.read_csv(filename)

    x_ticks = ("(" + df["num_nodes"].astype("str") + ", " + df["num_thread"].astype("str")) + ")"

    plt.plot(x_ticks, df["execution-time"].iloc[0] / df["execution-time"], label=get_arguments_from_filename(filename))

plt.legend()
plt.tight_layout()
plt.savefig("speedup.png")
