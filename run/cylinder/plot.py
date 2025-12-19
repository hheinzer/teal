from re import U
import subprocess
from glob import glob

import pyvista as pv
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable


def plot(fname):
    mesh = pv.read(fname)
    mesh = mesh.slice(normal="z")

    scalars = "velocity"
    cmap = "RdYlBu_r"
    vmin, vmax = 0.0, 1.7

    plot = pv.Plotter(lighting="three lights", off_screen=True)
    plot.add_mesh(mesh, scalars=scalars, cmap=cmap, clim=(vmin, vmax), show_scalar_bar=False)
    plot.camera.tight(view="xy")
    img = plot.screenshot(transparent_background=True)
    plot.close()

    plt.style.use("dark_background")
    fig, ax = plt.subplots()

    im = ax.imshow(img, extent=mesh.bounds, vmin=vmin, vmax=vmax, cmap=cmap)
    ax.set_xlabel("x-axis")
    ax.set_ylabel("y-axis")

    divider = make_axes_locatable(ax)
    cax = divider.append_axes("top", size=0.1, pad=0.1)
    cbar = plt.colorbar(im, cax=cax, orientation="horizontal")
    cbar.set_label(scalars)

    cax.xaxis.set_ticks_position("top")
    cax.xaxis.set_label_position("top")
    cax.tick_params(axis="x", top=True, labeltop=True, bottom=False, labelbottom=False)

    fig.tight_layout()
    plt.savefig(f"{fname}.png", dpi=200, bbox_inches="tight")
    plt.close()


def main():
    prefix = "bin/cylinder/run"
    fnames = sorted(glob(f"{prefix}_0*.vtkhdf"))
    for fname in fnames:
        plot(fname)

    subprocess.run(
        [
            "ffmpeg",
            "-y",
            "-framerate",
            "20",
            "-i",
            f"{prefix}_%05d.vtkhdf.png",
            "-pix_fmt",
            "yuv420p",
            "-vf",
            "pad=ceil(iw/2)*2:ceil(ih/2)*2",
            "-r",
            "60",
            f"{prefix}.mp4",
        ]
    )


if __name__ == "__main__":
    main()
