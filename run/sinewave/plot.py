import subprocess
from glob import glob

import numpy as np
import pyvista as pv


def plot(fname, azimuth):
    mesh = pv.read(fname)
    mesh = mesh.threshold((0, 0), scalars="entity")

    scalars = "density"
    cmap = "RdYlBu_r"
    vmin, vmax = 1.9, 2.1

    pv.set_plot_theme("dark")
    plot = pv.Plotter(lighting="three lights", off_screen=True, window_size=(1024, 1024))
    plot.add_mesh(mesh, scalars=scalars, cmap=cmap, clim=(vmin, vmax))

    plot.camera_position = "xy"
    plot.camera.azimuth = azimuth
    plot.camera.elevation += 30
    plot.camera.focal_point = np.array(plot.camera.focal_point) - (0, 0.3, 0)
    plot.screenshot(f"{fname}.png")


def main():
    prefix = "bin/sinewave/navier_stokes"
    fnames = sorted(glob(f"{prefix}_0*.vtkhdf"))
    azimuth = 30
    for fname in fnames:
        plot(fname, azimuth)
        azimuth += 45 / len(fnames)

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
