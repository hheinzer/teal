from glob import glob

import numpy as np
import pyvista as pv


def main():
    fname = sorted(glob("bin/transonic_wing/run_0*.vtkhdf"))[-1]
    mesh = pv.read(fname)

    wing = mesh.threshold((1, 1), scalars="entity")

    domain = mesh.threshold((0, 0), scalars="entity")
    domain = domain.cell_data_to_point_data()
    domain = domain.compute_derivative(scalars="velocity", gradient=False, qcriterion=True)
    domain["mach"] = np.linalg.norm(domain["velocity"], axis=1) / np.sqrt(
        mesh["heat capacity ratio"] * domain["pressure"] / domain["density"]
    )

    contour = domain.contour([5], scalars="qcriterion", compute_normals=True)

    streamlines = domain.streamlines(
        vectors="velocity", source_center=(0, 0, 1.4), source_radius=0.25, n_points=15
    )
    tubes = streamlines.tube(radius=0.005)
    tubes = tubes.clip_box(bounds=(0, 2, -1, 1, 0, 2), invert=False)

    plot = pv.Plotter(lighting="three lights")
    plot.add_mesh(wing, scalars="pressure", cmap="Blues", smooth_shading=True)
    plot.add_mesh(contour, scalars="mach", cmap="Reds", smooth_shading=True, opacity=0.75)
    plot.add_mesh(tubes, scalars="mach", cmap="Reds", smooth_shading=True)

    plot.show_axes()
    plot.camera_position = "xy"
    plot.camera.azimuth -= 30
    plot.camera.elevation += 30
    plot.camera.focal_point = np.array(plot.camera.focal_point) - (0, 0.3, 0)
    plot.camera.zoom(1.3)
    plot.save_graphic("bin/transonic_wing/transonic_wing.pdf")
    plot.show()


if __name__ == "__main__":
    main()
