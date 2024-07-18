from argparse import ArgumentParser
from os.path import splitext, basename
import numpy as np
import matplotlib.pyplot as plt

parser = ArgumentParser()
parser.add_argument("fnames", nargs="+", help="files that contain benchmark")
args = parser.parse_args()

fig, ax1 = plt.subplots()
ax2 = ax1.twinx()

max_ncores = 0
for fname in args.fnames:
    label = splitext(basename(fname))[0]
    label = label.replace("_", " ").replace("benchmark", "")
    label = " ".join(label.split())

    data = np.genfromtxt(fname)
    data = data[data[:, 0].argsort()]
    ncores, wtime = data[:, 0], data[:, 1]
    max_ncores = max(max_ncores, ncores[-1])

    # https://scicomp.ethz.ch/wiki/Parallel_efficiency
    ax1.plot(ncores, wtime[0] / wtime, label=label)
    ax2.plot(ncores, wtime[0] / (ncores * wtime))

ax1.plot([1, max_ncores], [1, max_ncores], "--", label="ideal")
ax2.plot([1, max_ncores], [1, 1], "--")

ax2.set_ylim(-0.1, 1.1)

ax1.legend(loc="lower right")
ax1.set_xlabel("number of processors")
ax1.set_ylabel("speedup")
ax2.set_ylabel("efficiency")

fig.tight_layout()
plt.show()
