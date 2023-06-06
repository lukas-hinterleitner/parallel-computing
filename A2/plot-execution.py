import pandas as pd
import numpy as np
import glob
import matplotlib.pyplot as plt


def get_arguments_from_filename(f: str) -> str:
    array = f.replace(".csv", "").split("-")

    i = 1 if ("omp" in f) else 0

    n = array[4 + i]
    m = array[3 + i]
    max_iterations = array[5 + i]
    epsilon = array[6 + i]

    return f"m={m}, n={n}, max_iter={max_iterations}, \u03B5={epsilon}"

omp_static_execution_time_files = glob.glob("execution-time-omp-static-*")
omp_dynamic_execution_time_files = glob.glob("execution-time-omp-dynamic-*")
mpi_execution_time_files = glob.glob("execution-time-mpi-*")

plt.figure(figsize=(15, 4))

plt.grid()
plt.title("OMP speedup")
plt.xlabel("number of threads")
plt.ylabel("speedup")
plt.xticks(np.arange(1, 33))
plt.yticks(np.arange(1, 20))

for filename in sorted(omp_static_execution_time_files):
    df = pd.read_csv(filename)
    plt.plot(df["num_thread"], df["execution-time"].iloc[0] / df["execution-time"], label=get_arguments_from_filename(filename))

plt.legend()
plt.tight_layout()
plt.savefig("omp-static-speedup.png")
plt.clf()

plt.grid()
plt.title("OMP speedup")
plt.xlabel("number of threads")
plt.ylabel("speedup")
plt.xticks(np.arange(1, 33))
plt.yticks(np.arange(1, 20))

for filename in sorted(omp_dynamic_execution_time_files):
    df = pd.read_csv(filename)
    plt.plot(df["num_thread"], df["execution-time"].iloc[0] / df["execution-time"], label=get_arguments_from_filename(filename))

plt.legend()
plt.tight_layout()
plt.savefig("omp-dynamic-speedup.png")
plt.clf()

plt.title("MPI speedup")
plt.xlabel("configuration (nodes, processes)")
plt.ylabel("speedup")
plt.grid()
for filename in sorted(mpi_execution_time_files):
    df = pd.read_csv(filename)

    x_ticks = ("(" + df["num_nodes"].astype("str") + ", " + df["num_thread"].astype("str")) + ")"

    plt.plot(x_ticks, df["execution-time"].iloc[0] / df["execution-time"], label=get_arguments_from_filename(filename))

plt.legend()
plt.tight_layout()
plt.savefig("mpi-speedup.png")
plt.clf()
