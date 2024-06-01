import numpy as np
import matplotlib.pyplot as plt
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument("fnames", nargs="+", help="files that contains the benchmark output")
args = parser.parse_args()

fig, ax = plt.subplots()

max_ncores = 0
for fname in args.fnames:
    label = "".join(fname.split("/")[-1].split(".")[:-1])
    label = label.replace("benchmark_", "").replace("_", " ")
    data = np.genfromtxt(fname)
    data = data[data[:, 0].argsort()]
    ncores = data[:, 0]
    wtime = data[:, 1]
    max_ncores = max(max_ncores, ncores[-1])
    ax.plot(ncores, wtime[0] / wtime, label=label)

ax.plot([1, max_ncores], [1, max_ncores], "--", label="ideal")

ax.set_xlabel("number of processors")
ax.set_ylabel("speedup")
ax.legend()
ax.set_aspect(1.0 / ax.get_data_ratio())

fig.tight_layout()
fig.savefig("bin/benchmark.pdf", bbox_inches="tight")
