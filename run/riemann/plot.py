from concurrent.futures import ProcessPoolExecutor
import os
import os.path
import subprocess

import matplotlib.pyplot as plt
from pypdf import PdfWriter
import pyvista as pv

CASES = ["sod", "toro1", "toro2", "toro3", "toro4", "toro5"]
FLUXES = ["godunov", "roe", "hll", "hllc", "hlle", "ausmd", "ausmdv", "lxf"]
ORDERS = ["1", "2"]
FIELDS = ["density", "velocity", "pressure", "energy"]


def ensure_solution(bin, flux, order):
    fname = f"{bin}_{flux}_{order}_00001.hdf"
    if not os.path.isfile(fname):
        subprocess.run([bin, "-q", flux, order], check=True)
    return fname


def read_solution(fname):
    mesh = pv.read(fname)
    mesh.cell_data.remove("vtkGhostType")
    gamma = mesh.field_data["heat capacity ratio"][0]

    def energy(density, velocity, pressure):
        return pressure / (gamma - 1) + 0.5 * density * velocity[:, 0] ** 2

    mesh["energy"] = energy(mesh["density"], mesh["velocity"], mesh["pressure"])
    mesh["energy-ref"] = energy(mesh["density-ref"], mesh["velocity-ref"], mesh["pressure-ref"])

    mesh = mesh.sample_over_line((0, 0.5, 0.5), (1, 0.5, 0.5), resolution=1000)
    idx = mesh["Distance"].argsort()

    def column(name):
        data = mesh[name]
        return data[idx] if data.ndim == 1 else data[:, 0][idx]

    distance = mesh["Distance"][idx]
    values = {name: column(name) for name in FIELDS + [f"{f}-ref" for f in FIELDS]}
    return distance, values


def plot(args):
    bin, order = args
    fig, axs = plt.subplots(nrows=2, ncols=2, layout="constrained")

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

    fname = f"{bin}_{order}.pdf"
    plt.savefig(fname)
    return fname


def main():
    subprocess.run(["make"], check=True)
    tasks = [(f"bin/riemann/{case}", order) for case in CASES for order in ORDERS]

    with ProcessPoolExecutor() as pool:
        fnames = list(pool.map(plot, tasks))

    writer = PdfWriter()
    for fname in fnames:
        writer.append(fname)
    with open("bin/riemann/riemann.pdf", "wb") as f:
        writer.write(f)

    for fname in fnames:
        os.remove(fname)


if __name__ == "__main__":
    main()
