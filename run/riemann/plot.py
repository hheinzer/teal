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
    fname = f"{bin}_{flux}_{order}_00001.vtkhdf"
    if not os.path.isfile(fname):
        subprocess.run([bin, "-q", flux, order], check=True)
    return fname


def read_solution(fname):
    mesh = pv.read(fname)
    mesh = mesh.sample_over_line((0, 0.5, 0.5), (1, 0.5, 0.5), resolution=1000)
    idx = mesh["Distance"].argsort()

    def column(name):
        data = mesh[name]
        return data[idx] if data.ndim == 1 else data[:, 0][idx]

    distance = mesh["Distance"][idx]
    values = {name: column(name) for name in FIELDS + [f"exact {f}" for f in FIELDS]}
    return distance, values


def plot_order(bin, order):
    fig, axs = plt.subplots(nrows=2, ncols=2)

    for i, flux in enumerate(FLUXES):
        fname = ensure_solution(bin, flux, order)
        distance, data = read_solution(fname)

        for ax, field in zip(axs.flat, FIELDS):
            if i == 0:
                ax.plot(distance, data[f"exact {field}"], "k", label="exact")

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
