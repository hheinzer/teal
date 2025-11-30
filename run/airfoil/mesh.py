import argparse
from pathlib import Path

import gmsh
import numpy as np


def parse_args():
    parser = argparse.ArgumentParser(description="Create an airfoil mesh")
    parser.add_argument("fname", type=Path, help="airfoil coordinate file (x y)")
    parser.add_argument("--upper", type=int, default=200, help="points along the upper surface")
    parser.add_argument("--lower", type=int, default=200, help="points along the lower surface")
    parser.add_argument("--outer", type=int, default=100, help="points on the farfield circle")
    parser.add_argument("--radius", type=float, default=20.0, help="farfield circle radius")
    parser.add_argument("--gui", action="store_true", help="launch the Gmsh GUI")
    return parser.parse_args()


def read_points(fname):
    points = np.genfromtxt(fname, skip_header=1, dtype=float)
    points = np.atleast_2d(points)

    if points.shape[1] > 2:
        points = points[:, :2]
    elif points.shape[1] < 2:
        raise ValueError("Need at least 2 columns (x y)")

    if np.allclose(points[0], points[-1]):
        points = points[:-1]

    return points


def get_lateral_extent(tag):
    xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.getBoundingBox(2, tag)
    if abs(zmax - zmin) < 1e-12:
        return -1.0  # cap, not lateral
    return max(abs(xmin), abs(xmax), abs(ymin), abs(ymax))


def create_model(points, upper, lower, outer, radius):
    gmsh.clear()
    gmsh.model.add("airfoil")

    for tag, (x, y) in enumerate(points, start=1):
        gmsh.model.occ.addPoint(x, y, 0.0, tag=tag)

    tag_back = np.argmax(points[:, 0]) + 1
    tag_front = np.argmin(points[:, 0]) + 1
    tag_last = points.shape[0]

    if not (1 <= tag_back <= tag_front <= tag_last):
        raise ValueError("Point ordering does not match expected TE->upper->LE->lower->TE layout")

    tag_upper = gmsh.model.occ.addSpline(list(range(tag_back, tag_front + 1)))
    tag_lower = gmsh.model.occ.addSpline(list(range(tag_front, tag_last + 1)) + [tag_back])

    tag_outer = gmsh.model.occ.addCircle(0.5, 0.0, 0.0, radius)

    curve_outer = gmsh.model.occ.addCurveLoop([tag_outer])
    curve_inner = gmsh.model.occ.addCurveLoop([tag_upper, tag_lower])
    plane = gmsh.model.occ.addPlaneSurface([curve_outer, curve_inner])

    extrude = gmsh.model.occ.extrude(
        [(2, plane)], 0.0, 0.0, 1.0, numElements=[1], heights=[1.0], recombine=True
    )
    domain = [tag for (dim, tag) in extrude if dim == 3]

    if not domain:
        raise RuntimeError("Extrude did not create a volume")

    gmsh.model.occ.synchronize()

    gmsh.model.mesh.setTransfiniteCurve(tag_upper, upper + 1)
    gmsh.model.mesh.setTransfiniteCurve(tag_lower, lower + 1)
    gmsh.model.mesh.setTransfiniteCurve(tag_outer, outer + 1)

    boundaries = gmsh.model.getBoundary([(3, domain[0])])
    surfaces = [tag for (dim, tag) in boundaries if dim == 2]

    extent = [tag for tag in surfaces if get_lateral_extent(tag) >= 0]
    farfield = max(extent, key=get_lateral_extent)
    airfoil = [surface for surface in extent if surface != farfield]

    gmsh.model.addPhysicalGroup(2, airfoil, name="airfoil")
    gmsh.model.addPhysicalGroup(2, [farfield], name="farfield")
    gmsh.model.addPhysicalGroup(3, domain, name="domain")


def main():
    args = parse_args()
    points = read_points(args.fname)

    gmsh.initialize()
    create_model(points, args.upper, args.lower, args.outer, args.radius)
    gmsh.model.mesh.generate(3)

    if args.gui:
        gmsh.fltk.run()

    gmsh.write(str(args.fname.with_suffix(".msh")))
    gmsh.finalize()


if __name__ == "__main__":
    main()
