from argparse import ArgumentParser, BooleanOptionalAction as Boolean
import h5py
import numpy as np
import pyvista as pv

parser = ArgumentParser(description="Visualize an unstructured grid")
parser.add_argument("fname", help="Input file")
parser.add_argument("field", nargs="?", default=None, help="Field to plot")
parser.add_argument("-l", "--list", action=Boolean, help="List all plottable fields")
parser.add_argument("-p", "--points", action=Boolean, help="Interpolate to points")
parser.add_argument("-x", "--explode", action=Boolean, help="Explode cells")
parser.add_argument("-c", "--clip", action=Boolean, help="Add clip plane")
parser.add_argument("-e", "--edges", action=Boolean, help="Show edges")
parser.add_argument("--crinkle", action=Boolean, help="Crinkle clip")
parser.add_argument("--cmap", default="coolwarm", help="Colormap")
args = parser.parse_args()


def read_mesh(fname):
    with h5py.File(fname, "r") as file:
        coords = file["nodes/coord"][...]
        node_offs = file["cells/node/off"][...]
        node_idxs = file["cells/node/idx"][...]
        types = file["cells/type"][...]
        fields = {
            name: field[...]
            for name, field in file["cells/"].items()
            if isinstance(field, h5py.Dataset) and field.shape[0] == types.shape[0]
        }

    cells = []
    for beg, end in zip(node_offs, node_offs[1:]):
        cells.append(end - beg)
        cells.extend(node_idxs[beg:end])

    mesh = pv.UnstructuredGrid(cells, types, coords)
    for name, field in fields.items():
        mesh[name] = field

    return mesh


def main():
    mesh = read_mesh(args.fname) if "_mesh.h5" in args.fname else pv.read(args.fname)
    mesh.set_active_scalars(args.field)

    if args.list:
        print("\n".join(mesh.array_names))
        return

    if args.points:
        mesh = mesh.cell_data_to_point_data()

    pvp = pv.Plotter(lighting="three lights")

    if args.explode:
        mesh = mesh.explode()
        pvp.disable_ssao()
    else:
        vol_types = np.array([10, 12, 13, 14], dtype=np.uint8)
        keep = np.where(np.isin(mesh.celltypes, vol_types))[0]
        mesh = mesh.extract_cells(keep).clean()

    if args.clip:
        pvp.add_mesh_clip_plane(mesh, show_edges=args.edges, crinkle=args.crinkle, cmap=args.cmap)
    else:
        pvp.add_mesh(mesh, show_edges=args.edges, cmap=args.cmap)

    pvp.show_bounds()
    pvp.view_xy()
    pvp.show()


if __name__ == "__main__":
    main()
