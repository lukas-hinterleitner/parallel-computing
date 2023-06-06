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
plt.title("OMP speedup with static iteration distribution")
plt.xlabel("number of threads")
plt.ylabel("speedup")
plt.xticks(np.arange(1, 33))
plt.yticks(np.arange(1, 20))

for filename in sorted(omp_static_execution_time_files):
    df = pd.read_csv(filename)
    df["speedup"] = df["execution-time"].iloc[0] / df["execution-time"]
    plt.plot(df["num_thread"], df["speedup"], label=get_arguments_from_filename(filename))

    df.to_csv(filename, index=False)

plt.legend()
plt.tight_layout()
plt.savefig("omp-static-speedup.png")
plt.clf()

plt.grid()
plt.title("OMP speedup with dynamic iteration distribution")
plt.xlabel("number of threads")
plt.ylabel("speedup")
plt.xticks(np.arange(1, 33))
plt.yticks(np.arange(1, 20))

for filename in sorted(omp_dynamic_execution_time_files):
    df = pd.read_csv(filename)
    df["speedup"] = df["execution-time"].iloc[0] / df["execution-time"]
    plt.plot(df["num_thread"], df["speedup"], label=get_arguments_from_filename(filename))

    df.to_csv(filename, index=False)

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
    df["speedup"] = df["execution-time"].iloc[0] / df["execution-time"]
    plt.plot(x_ticks, df["speedup"], label=get_arguments_from_filename(filename))
    df.to_csv(filename, index=False)

plt.legend()
plt.tight_layout()
plt.savefig("mpi-speedup.png")
plt.clf()


plt.title("MPI/OMP single node speedup")
plt.xlabel("processes")
plt.ylabel("speedup")
plt.grid()
for mpi_filename, omp_filename in zip(sorted(mpi_execution_time_files), sorted(omp_static_execution_time_files)):
    mpi_df = pd.read_csv(mpi_filename)
    omp_df = pd.read_csv(omp_filename)

    mpi_df["speedup"] = mpi_df["execution-time"].iloc[0] / mpi_df["execution-time"]
    omp_df["speedup"] = omp_df["execution-time"].iloc[0] / omp_df["execution-time"]

    filtered_mpi_df = mpi_df[mpi_df["num_nodes"] == 1]
    filtered_omp_df = omp_df[omp_df["num_thread"].isin(filtered_mpi_df["num_thread"])]

    plt.plot(filtered_mpi_df["num_thread"], filtered_mpi_df["speedup"], label=f"MPI: {get_arguments_from_filename(mpi_filename)}")
    plt.plot(filtered_omp_df["num_thread"], filtered_omp_df["speedup"], label=f"OMP: {get_arguments_from_filename(omp_filename)}")

plt.legend()
plt.tight_layout()
plt.savefig("mpi-omp-single-node-speedup.png")
plt.clf()
