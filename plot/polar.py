import h5py
import numpy as np
import matplotlib.pyplot as plt
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument("fname", help="h5 file that contains the polar")
parser.add_argument("refs", nargs="*", help="dat files that contain reference polars")
args = parser.parse_args()

fig, ax = plt.subplots(nrows=2, height_ratios=(3, 1))

data = h5py.File(args.fname, "r")
x = data["x"][...]
cp = data["cp"][...]
cl = data["cl"][0]
cd = data["cd"][0]

mean = np.mean(x, axis=0)
theta = np.arctan2(x[:, 1] - mean[1], x[:, 0] - mean[0])
order = np.argsort(theta)
order = np.concatenate((order, [order[0]]))

ax[0].plot(x[order, 0], cp[order])
ax[1].plot(x[order, 0], x[order, 1])

if args.refs:
    for ref in args.refs:
        data = np.genfromtxt(ref, skip_header=1)
        x = data[:, 0]
        cp = data[:, 1]
        ax[0].plot(x, cp, "--", label=ref.split("_")[-1].replace(".dat", ""))
    ax[0].legend(title="reference")

ax[0].invert_yaxis()
ax[0].set_xlabel("$x$")
ax[0].set_ylabel("$c_p$")
ax[0].set_title(f"$c_l$ = {cl:.4f}, $c_d$ = {cd:.4f}")
ax[1].set_aspect("equal")
ax[1].set_axis_off()

fig.tight_layout()
fig.savefig(args.fname.replace(".h5", ".pdf"), bbox_inches="tight")
