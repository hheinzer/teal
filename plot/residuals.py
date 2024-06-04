import numpy as np
import matplotlib.pyplot as plt
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument("fname", help="file that contains the residuals")
args = parser.parse_args()

fig, ax = plt.subplots()

with open(args.fname) as file:
    name = file.readline().split()

data = np.genfromtxt(args.fname, skip_header=1)
for v in range(data.shape[1]):
    if all(np.isclose(data[:, v], 0)):
        continue
    ax.plot(data[:, v], label=name[v])

ax.set_xlabel("iteration")
ax.set_ylabel("residual")
ax.legend()
ax.set_yscale("log")

fig.tight_layout()
fig.savefig(args.fname.replace(".dat", ".pdf"), bbox_inches="tight")
