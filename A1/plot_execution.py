import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

sequential_execution_time = 10.1871

plt.figure(figsize=(15,5))
plt.xlabel("number of threads")
plt.ylabel("speedup")
plt.grid(axis="y")
plt.title("execution on ALMA cluster of University of Vienna")

for method in ["mutex", "atomic"]:
    filename = f"execution_time_{method}.csv"

    execution_times = pd.read_csv(filename).to_numpy()

    num_threads = execution_times[:, 0]
    speedup = sequential_execution_time / execution_times[:, 1]

    plt.xticks(np.arange(0, max(num_threads)+1))
    plt.plot(num_threads, speedup, label=method)
    plt.legend()

plt.savefig("speedup.png")
