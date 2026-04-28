import os.path
import subprocess
from concurrent.futures import ProcessPoolExecutor

import matplotlib.pyplot as plt
import pyvista as pv

CASES = ["sod", "toro1", "toro2", "toro3", "toro4", "toro5"]
ORDERS = ["1", "2"]
FLUXES = ["godunov", "roe", "hll", "hllc", "hlle", "ausmd", "ausmdv", "lxf"]
FIELDS = ["density", "velocity", "pressure", "energy"]


def ensure_solution(bin, flux, order):
    fname = f"{bin}_{flux}_{order}_00001.hdf"
    if not os.path.isfile(fname):
        subprocess.run([bin, "-q", flux, order], check=True)
    return fname


def read_solution(fname):
    mesh = pv.read(fname)
    gamma = mesh.field_data["heat capacity ratio"][0]

    def energy(density, velocity, pressure):
        return pressure / (gamma - 1) + 0.5 * density * velocity[:, 0] ** 2

    mesh["energy"] = energy(mesh["density"], mesh["velocity"], mesh["pressure"])
    mesh["energy-ref"] = energy(mesh["density-ref"], mesh["velocity-ref"], mesh["pressure-ref"])

    def column(name):
        data = mesh[name]
        return data if data.ndim == 1 else data[:, 0]

    centers = mesh["center"][:, 0]
    idx = centers.argsort()
    distance = centers[idx]
    values = {name: column(name)[idx] for name in FIELDS + [f"{f}-ref" for f in FIELDS]}
    return distance, values


def plot_order(bin, order):
    fig, axs = plt.subplots(nrows=2, ncols=2)

    for i, flux in enumerate(FLUXES):
        fname = ensure_solution(bin, flux, order)
        distance, data = read_solution(fname)

        for ax, field in zip(axs.flat, FIELDS):
            if i == 0:
                ax.plot(distance, data[f"{field}-ref"], "k", label="exact")

            ax.plot(distance, data[field], label=f"{flux}")
            ax.set_xlabel("x-axis")
            ax.set_ylabel(field)

    axs[0, 1].legend(fontsize="small", loc="upper left", bbox_to_anchor=(1, 1))

    fig.suptitle(f"{bin.split("/")[-1]} O{order}")
    fig.tight_layout()
    plt.savefig(f"{bin}_{order}.pdf")


def run_plot(args):
    bin, order = args
    plot_order(bin, order)
    return bin, order


def main():
    subprocess.run(["make"], check=True)
    tasks = [(f"bin/riemann/{case}", order) for case in CASES for order in ORDERS]
    with ProcessPoolExecutor() as pool:
        for bin, order in pool.map(run_plot, tasks):
            print(bin, order)


if __name__ == "__main__":
    main()
