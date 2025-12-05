from pathlib import Path

import numpy as np
import pyvista as pv
import matplotlib.pyplot as plt


def read_mesh(fname):
    mesh = pv.read(fname)
    if isinstance(mesh, pv.PartitionedDataSet):
        blocks = pv.MultiBlock([part for part in mesh if part is not None])
        mesh = blocks.combine()
    return mesh


def sample(fname, pointa, pointb):
    name = " ".join(Path(fname).stem.split("_")[:-1])
    mesh = read_mesh(fname)
    mesh = mesh.sample_over_line(pointa, pointb, resolution=1000)
    return name, mesh.point_data


def main():
    datas = [
        sample("bin/explosion/1D_alpha_0_00001.vtkhdf", (0, 0.5, 0.5), (1, 0.5, 0.5)),
        sample("bin/explosion/1D_alpha_1_00001.vtkhdf", (0, 0.5, 0.5), (1, 0.5, 0.5)),
        sample("bin/explosion/1D_alpha_2_00001.vtkhdf", (0, 0.5, 0.5), (1, 0.5, 0.5)),
        sample("bin/explosion/2D_00001.vtkhdf", (0, 0, 0), (1, 1, 0)),
        sample("bin/explosion/3D_00001.vtkhdf", (0, 0, 0), (1, 1, 1)),
    ]

    fig, axs = plt.subplots(nrows=2, ncols=2)

    keys = ["density", "velocity", "pressure", "energy"]
    for ax, key in zip(axs.flat, keys):
        for data in datas:
            label = data[0]
            dist = data[1]["Distance"]
            data = data[1][key]
            if data.ndim == 2:
                data = np.linalg.norm(data, axis=1)
            ax.plot(dist[dist <= 1], data[dist <= 1], label=label)
            ax.set_xlabel("distance")
            ax.set_ylabel(key)

    axs[0, 0].legend(fontsize="small")

    fig.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()
