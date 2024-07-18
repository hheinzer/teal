from argparse import ArgumentParser
from os.path import basename
import pyvista as pv
import matplotlib.pyplot as plt

parser = ArgumentParser()
parser.add_argument("fnames", nargs="+", help="files that should be plotted")
parser.add_argument("scalar", help="scalar to plot")
parser.add_argument("-a", "--pointa", help="x coordinate of start point", type=float)
parser.add_argument("-b", "--pointb", help="x coordinate of end point", type=float)
parser.add_argument("-r", "--resolution", help="line resolution", type=int)
args = parser.parse_args()

fig, ax = plt.subplots()

for fname in args.fnames:
    mesh = pv.read(fname)

    xmin = mesh.points.min(axis=0)
    xmax = mesh.points.max(axis=0)
    pointa = 0.5 * (xmin + xmax)
    pointb = 0.5 * (xmin + xmax)
    pointa[0] = args.pointa if args.pointa is not None else xmin[0]
    pointb[0] = args.pointb if args.pointb is not None else xmax[0]
    resolution = args.resolution if args.resolution is not None else mesh.n_cells + 1
    sample = mesh.sample_over_line(pointa, pointb, resolution)

    fname = " ".join(basename(fname).split("_")[:-1])
    for scalar in sample.array_names:
        if args.scalar in scalar:
            label = (scalar + fname).replace(args.scalar, "")
            zorder = 1 if "exact" in scalar else 2
            ax.plot(sample["Distance"], sample[scalar], label=label, zorder=zorder)

ax.legend(loc="upper left", bbox_to_anchor=(1, 1))
ax.set_xlabel("distance")
ax.set_ylabel(args.scalar)

fig.tight_layout()
plt.show()
