from glob import glob

import pyvista as pv
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable


def main():
    fname = sorted(glob("bin/forward_facing_step/run_0*.vtkhdf"))[-1]
    mesh = pv.read(fname)
    mesh = mesh.slice(normal="z")

    scalars = "density"
    cmap = "RdYlBu_r"

    plot = pv.Plotter(lighting="three lights", off_screen=True)
    plot.add_mesh(mesh, scalars=scalars, cmap=cmap, show_scalar_bar=False)
    plot.camera.tight(view="xy")
    img = plot.screenshot()
    plot.close()

    fig, ax = plt.subplots()

    vmin, vmax = mesh[scalars].min(), mesh[scalars].max()
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
    plt.savefig("bin/forward_facing_step/forward_facing_step.pdf", bbox_inches="tight")
    plt.show()


if __name__ == "__main__":
    main()
