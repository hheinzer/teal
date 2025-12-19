import subprocess
from glob import glob

import pyvista as pv
from mpi4py import MPI


def plot(fname, azimuth):
    mesh = pv.read(fname)
    mesh = mesh.threshold((0, 0), scalars="entity")
    corners = mesh.outline_corners()

    mesh = mesh.cell_data_to_point_data()
    mesh = mesh.compute_derivative(scalars="velocity", gradient=False, qcriterion=True)
    mesh = mesh.contour([0.3], scalars="qcriterion", compute_normals=True)

    scalars = "velocity"
    cmap = "RdYlBu_r"
    vmin, vmax = 0.0, 1.0

    pv.set_plot_theme("dark")
    plot = pv.Plotter(off_screen=True, window_size=(1028, 1028))
    plot.add_mesh(corners)
    plot.add_mesh(mesh, scalars=scalars, cmap=cmap, clim=(vmin, vmax), smooth_shading=True)

    plot.camera_position = "xy"
    plot.camera.azimuth = azimuth
    plot.camera.elevation = 30
    plot.screenshot(f"{fname}.png")
    plot.close()


def main():
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()

    prefix = "bin/taylor_green_vortex/run"
    fnames = sorted(glob(f"{prefix}_0*.vtkhdf"))

    for fname in fnames[rank::size]:
        azimuth = 30.0 + 45.0 * fnames.index(fname) / (len(fnames) - 1)
        plot(fname, azimuth)

    comm.Barrier()
    if rank == 0:
        subprocess.run(
            [
                "ffmpeg",
                "-y",
                "-framerate",
                "30",
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
