#!/bin/python
from argparse import ArgumentParser
from os.path import basename, commonprefix
import pyvista as pv
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np

parser = ArgumentParser()
parser.add_argument("fnames", nargs="+", help="files that should be plotted")
parser.add_argument("scalar", help="scalar to plot")
parser.add_argument("-x", help="x coordinates of slices")
parser.add_argument("-r", "--resolution", help="line resolution", type=int, default=100)
args = parser.parse_args()

fig, ax = plt.subplots()

linestyles = ["solid", "dashed", "dotted", "dashdot"]
for fname, ls in zip(args.fnames, linestyles):
    mesh = pv.read(fname)
    mesh = mesh.cell_data_to_point_data()
    xmin = np.array(mesh.bounds[0::2])
    xmax = np.array(mesh.bounds[1::2])
    delx = xmax - xmin

    X = (
        [float(x) for x in args.x.split()]
        if args.x is not None
        else np.arange(xmin[0], xmax[0], delx[0] / 5)
    )

    colors = mcolors.TABLEAU_COLORS
    for x, c in zip(X, colors):
        mean = None
        for z in np.arange(xmin[2], xmax[2], delx[2] / args.resolution):
            pointa = np.array([x, xmin[1], z])
            pointb = np.array([x, xmax[1], z])
            sample = mesh.sample_over_line(pointa, pointb, args.resolution)
            mean = mean + sample[args.scalar] if mean is not None else sample[args.scalar]
        mean /= args.resolution
        mean = np.linalg.norm(mean, axis=1) if len(mean.shape) > 1 else mean

        ax.plot(mean, sample.points[:, 1], c=c, ls=ls)

if len(args.fnames) > 1:
    for fname, ls in zip(args.fnames, linestyles):
        fname = basename(fname).replace("_", " ").replace(".hdf", "")
        ax.plot([], [], c="k", ls=ls, label=fname)

for x, c in zip(X, colors):
    ax.plot([], [], c=c, ls="", marker="s", label=f"x = {x:g}")

ax.legend(loc="upper left", bbox_to_anchor=(1, 1))
ax.set_xlabel(args.scalar)
ax.set_ylabel("y")

fig.tight_layout()
prefix = commonprefix([fname.replace(".hdf", "") for fname in args.fnames])
fig.savefig((prefix if prefix[-1] != "_" else prefix[:-1]) + "_channel.pdf", bbox_inches="tight")
