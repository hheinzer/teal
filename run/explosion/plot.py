import os.path
import subprocess

import matplotlib.pyplot as plt
import numpy as np
import pyvista as pv

FIELDS = ["density", "velocity", "pressure", "energy"]
CASES = [
    {
        "label": f"1D alpha {alpha}",
        "exe": "bin/explosion/1D",
        "args": [alpha],
        "file": f"bin/explosion/1D_alpha_{alpha}_00001.vtkhdf",
        "line": ((0, 0.5, 0.5), (1, 0.5, 0.5)),
    }
    for alpha in ["0", "1", "2"]
] + [
    {
        "label": "2D",
        "exe": "bin/explosion/2D",
        "args": [],
        "file": "bin/explosion/2D_00001.vtkhdf",
        "line": ((0, 0, 0), (1, 1, 0)),
    },
    {
        "label": "3D",
        "exe": "bin/explosion/3D",
        "args": [],
        "file": "bin/explosion/3D_00001.vtkhdf",
        "line": ((0, 0, 0), (1, 1, 1)),
    },
]


def ensure_solution(case):
    if not os.path.isfile(case["file"]):
        subprocess.run(["mpirun", "-n", "6", case["exe"], "-q", *case["args"]], check=True)
    return case["file"]


def sample_solution(case):
    fname = ensure_solution(case)
    mesh = pv.read(fname)
    mesh = mesh.sample_over_line(*case["line"], resolution=1000)
    idx = mesh["Distance"].argsort()
    mask = (0 <= mesh["Distance"]) & (mesh["Distance"] <= 1)
    idx = idx[mask[idx]]

    def column(name):
        data = mesh[name]
        if data.ndim == 1:
            return data[idx]
        if data.shape[1] == 1:
            return data[:, 0][idx]
        return data[idx]

    distance = mesh["Distance"][idx]
    values = {name: column(name) for name in FIELDS}
    if values["velocity"].ndim == 2:
        values["velocity"] = np.linalg.norm(values["velocity"], axis=1)
    return case["label"], distance, values


def plot(series):
    fig, axs = plt.subplots(nrows=2, ncols=2)

    for ax, field in zip(axs.flat, FIELDS):
        for label, distance, values in series:
            ax.plot(distance, values[field], label=label)
        ax.set_xlabel("distance")
        ax.set_ylabel(field)

    handles, labels = axs.flat[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=len(labels))
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    plt.savefig("bin/explosion/explosion.pdf", bbox_inches="tight")
    plt.show()


def main():
    subprocess.run(["make"], check=True)
    series = [sample_solution(case) for case in CASES]
    plot(series)


if __name__ == "__main__":
    main()
