from argparse import ArgumentParser

import h5py
import numpy as np
import matplotlib.pyplot as plt


def parse_args():
    parser = ArgumentParser(description="Visualize an airfoil polar")
    parser.add_argument("fname", help="input file")
    parser.add_argument("refs", nargs="*", help="reference polars")
    return parser.parse_args()


def get_order(center):
    mean = np.mean(center, axis=0)
    theta = np.arctan2(center[:, 1] - mean[1], center[:, 0] - mean[0])
    order = np.argsort(theta)
    return np.concatenate((order, [order[0]]))


def main():
    args = parse_args()

    with h5py.File(args.fname, "r") as file:
        center = file["center"][...]
        pressure_c = file["pressure_c"][...]
        lift_c = file["lift_c"][0]
        drag_c = file["drag_c"][0]

    order = get_order(center)
    center, pressure_c = center[order], pressure_c[order]

    fig, ax = plt.subplots(nrows=2, height_ratios=(3, 1))

    ax[0].plot(center[:, 0], pressure_c)
    ax[1].plot(center[:, 0], center[:, 1])

    for ref in args.refs:
        data = np.genfromtxt(ref, skip_header=True)
        ax[0].plot(data[:, 0], data[:, 1], "--")

    ax[0].set_xlabel("x-axis")
    ax[0].set_ylabel("pressure coefficient")
    ax[0].set_title(f"lift/drag coefficient = {lift_c:.4f} / {drag_c:.4f}")

    ax[0].invert_yaxis()
    ax[1].set_aspect("equal")
    ax[1].set_axis_off()

    fig.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()
