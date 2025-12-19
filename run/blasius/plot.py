from glob import glob

import numpy as np
import pyvista as pv
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.integrate import solve_ivp


def blasius(eta_max, n=2000, fpp0=0.332057336215):
    def rhs(_, y):
        f, fp, fpp = y
        fppp = -0.5 * f * fpp
        return [fp, fpp, fppp]

    eta = np.linspace(0, eta_max, n)
    y0 = [0, 0, fpp0]
    sol = solve_ivp(rhs, (eta[0], eta[-1]), y0, t_eval=eta)
    return eta, sol.y[1]


def main():
    fname = sorted(glob("bin/blasius/run_0*.vtkhdf"))[-1]
    mesh = pv.read(fname)
    mesh = mesh.clip(normal="y", origin=(0, 0.2, 0))
    mesh = mesh.slice(normal="z")

    scalars = "velocity"
    cmap = "RdYlBu_r"

    plot = pv.Plotter(lighting="three lights", off_screen=True)
    plot.add_mesh(mesh, scalars=scalars, cmap=cmap, show_scalar_bar=False)
    plot.camera.tight(view="xy")
    img = plot.screenshot()
    plot.close()

    _, x, ymin, ymax, zmin, zmax = mesh.bounds
    z = (zmin + zmax) / 2
    sample = mesh.sample_over_line((x, ymin, z), (x, ymax, z), resolution=1000)
    u = sample[scalars][:-1, 0]
    y = sample.points[:-1, 1]

    rho_inf, u_inf = 1, 1
    nu = mesh["dynamic viscosity"] / rho_inf
    eta = (y - ymin) * np.sqrt(u_inf / (nu * x))
    eta_tab, fp_tab = blasius(max(eta.max(), 10))
    fp = np.interp(eta, eta_tab, fp_tab, left=0.0, right=1.0)
    u_blasius = u_inf * fp

    fig, ax = plt.subplots(nrows=2, height_ratios=[1, 2])

    mag = np.linalg.norm(mesh[scalars], axis=1)
    im = ax[0].imshow(img, extent=mesh.bounds, vmin=mag.min(), vmax=mag.max(), cmap=cmap)
    ax[0].set_xlabel("x-axis")
    ax[0].set_ylabel("y-axis")

    divider = make_axes_locatable(ax[0])
    cax = divider.append_axes("top", size=0.1, pad=0.1)
    cbar = plt.colorbar(im, cax=cax, orientation="horizontal")
    cbar.set_label("|velocity|")

    cax.xaxis.set_ticks_position("top")
    cax.xaxis.set_label_position("top")
    cax.tick_params(axis="x", top=True, labeltop=True, bottom=False, labelbottom=False)

    ax[1].plot(u, y, label="teal")
    ax[1].plot(u_blasius, y, "--", label="blasius")
    ax[1].set_xlabel("x-velocity")
    ax[1].set_ylabel("y-axis")
    ax[1].legend()

    fig.subplots_adjust(hspace=0.3)
    plt.savefig("bin/blasius/blasius.pdf", bbox_inches="tight")
    plt.show()


if __name__ == "__main__":
    main()
