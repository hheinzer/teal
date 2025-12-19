from glob import glob

import pyvista as pv
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable


def main():
    fname = sorted(glob("bin/lid_driven_cavity/run_0*.vtkhdf"))[-1]
    mesh = pv.read(fname)
    mesh = mesh.slice(normal="z")

    scalars = "velocity"
    cmap = "RdYlBu_r"
    vmin, vmax = 0.0, 1.0

    plot = pv.Plotter(lighting="three lights", off_screen=True)
    plot.add_mesh(mesh, scalars=scalars, cmap=cmap, clim=(vmin, vmax), show_scalar_bar=False)
    plot.camera.tight(view="xy")
    img = plot.screenshot()
    plot.close()

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
    plt.savefig("bin/lid_driven_cavity/lid_driven_cavity.pdf", bbox_inches="tight")
    plt.show()


if __name__ == "__main__":
    main()
