#!/bin/python
from argparse import ArgumentParser, BooleanOptionalAction
from os import remove, system
from os.path import commonprefix
from time import time

import matplotlib.pyplot as plt
from mpi4py import MPI
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
from pypdf import PdfWriter
import pyvista as pv

parser = ArgumentParser()
parser.add_argument("fnames", nargs="+", help="files that should be plotted")
parser.add_argument("scalar", help="scalar to plot")
parser.add_argument("--vmin", help="minimum scalar value", type=float)
parser.add_argument("--vmax", help="maximum scalar value", type=float)
parser.add_argument("-x", "--explode", help="plot exploded mesh", action=BooleanOptionalAction)
parser.add_argument("-e", "--edges", help="plot mesh edges", action=BooleanOptionalAction)
parser.add_argument("-p", "--points", help="plot over point data", action=BooleanOptionalAction)
parser.add_argument("-t", "--time", help="use time as title", action=BooleanOptionalAction)
parser.add_argument("-c", "--cmap", help="colormap", default="coolwarm")
args = parser.parse_args()


def main():
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()

    fnames = []
    for i in range(size):
        fnames += args.fnames[i::size]
    fnames = np.array_split(fnames, size)[rank]

    for fname in fnames:
        plot(fname)

    comm.Barrier()
    if rank == 0 and len(args.fnames) > 1:
        merge_pdfs(args.fnames, args.scalar)
        create_video(args.fnames, args.scalar)


def plot(fname):
    tstart = time()

    mesh = pv.read(fname)
    mesh = mesh.explode(1) if args.explode else mesh
    mesh = mesh.cell_data_to_point_data() if args.points else mesh
    extent = mesh.bounds
    vmin = (
        args.vmin
        if args.vmin is not None
        else (
            mesh[args.scalar].min()
            if len(mesh[args.scalar].shape) == 1
            else np.linalg.norm(mesh[args.scalar], axis=1).min()
        )
    )
    vmax = (
        args.vmax
        if args.vmax is not None
        else (
            mesh[args.scalar].max()
            if len(mesh[args.scalar].shape) == 1
            else np.linalg.norm(mesh[args.scalar], axis=1).max()
        )
    )

    pvp = pv.Plotter(lighting="three lights", off_screen=True, window_size=[4000, 4000])
    pvp.add_mesh(
        mesh,
        scalars=args.scalar,
        clim=(vmin, vmax),
        show_edges=args.edges,
        cmap=args.cmap,
        show_scalar_bar=False,
    )
    pvp.camera.tight()
    img = pvp.screenshot(None, transparent_background=True, return_img=True)
    pv.close_all()

    fig, ax = plt.subplots()

    im = ax.imshow(img, cmap=args.cmap, extent=extent, vmin=vmin, vmax=vmax, interpolation="none")

    div = make_axes_locatable(ax)
    cax = div.append_axes("right", size=0.15, pad=0.1)
    cb = plt.colorbar(im, cax=cax)

    ax.set_xlabel("x")
    ax.set_ylabel("y")
    cb.set_label(args.scalar)
    if args.time and "time" in mesh.array_names:
        ax.set_title(f"time = {mesh['time'][0]:16.4f}")

    fig.tight_layout()
    fig.savefig(fname.replace(".hdf", ".pdf"), bbox_inches="tight")
    fig.savefig(fname.replace(".hdf", ".png"), dpi=600, bbox_inches="tight")
    plt.close()

    tend = time()
    print(fname, tend - tstart, "s")


def merge_pdfs(fnames, scalar):
    with PdfWriter() as writer:
        for fname in fnames:
            writer.append(fname.replace(".hdf", ".pdf"))
        prefix = commonprefix(["_".join(fname.split("_")[:-1]) for fname in fnames])
        writer.write(f"{prefix}_{scalar}.pdf")
        writer.close()
    for fname in fnames:
        remove(fname.replace(".hdf", ".pdf"))


def create_video(fnames, scalar):
    prefix = commonprefix(["_".join(fname.split("_")[:-1]) for fname in fnames])
    system(
        " ".join(
            [
                "ffmpeg -y -loglevel warning",
                "-framerate 20",
                f"-i {prefix}_%05d.png",
                "-pix_fmt yuv420p",
                "-vf 'pad=ceil(iw/2)*2:ceil(ih/2)*2'",
                "-r 60",
                f"{prefix}_{scalar}.mp4",
            ]
        )
    )
    for fname in fnames:
        remove(fname.replace(".hdf", ".png"))


if __name__ == "__main__":
    main()
