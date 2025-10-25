from argparse import ArgumentParser, BooleanOptionalAction
import h5py
import numpy as np
import pyvista as pv

parser = ArgumentParser(description="Visualize an unstructured mesh")
parser.add_argument("fname", help="Input HDF5 file")
parser.add_argument("field", nargs="?", default=None, help="Field to plot")
parser.add_argument("-l", "--list", action=BooleanOptionalAction, help="List all plottable fields")
parser.add_argument("-x", "--explode", action=BooleanOptionalAction, help="Explode cells")
parser.add_argument("-e", "--edges", action=BooleanOptionalAction, help="Show edges")
parser.add_argument("-c", "--clip", action=BooleanOptionalAction, help="Add clip plane")
parser.add_argument("--crinkle", action=BooleanOptionalAction, help="Crinkle clip")
parser.add_argument("--cmap", default="coolwarm", help="Colormap")
args = parser.parse_args()

with h5py.File(args.fname, "r") as file:
    coords = file["nodes/coord"][...]
    node_offs = file["cells/node/off"][...]
    node_idxs = file["cells/node/idx"][...]
    types = file["cells/type"][...]

    if args.field:
        field = file[f"cells/{args.field}"][...]

    if args.list:
        for name, obj in file["cells/"].items():
            if isinstance(obj, h5py.Dataset):
                print(name)
        exit()

cells = []
for beg, end in zip(node_offs, node_offs[1:]):
    cells.append(end - beg)
    cells.extend(node_idxs[beg:end])
cells = np.array(cells, dtype=int)

mesh = pv.UnstructuredGrid(cells, types, coords)

if args.field:
    mesh[f"{args.field}"] = field

pvp = pv.Plotter(lighting="three lights")

if args.explode:
    mesh = mesh.explode()
    pvp.disable_ssao()
else:
    vol_types = np.array([10, 12, 13, 14], dtype=np.uint8)
    keep = np.where(np.isin(types, vol_types))[0]
    mesh = mesh.extract_cells(keep).clean()

if args.clip:
    pvp.add_mesh_clip_plane(mesh, show_edges=args.edges, cmap=args.cmap, crinkle=args.crinkle)
else:
    pvp.add_mesh(mesh, show_edges=args.edges, cmap=args.cmap)

pvp.show_bounds()
pvp.view_xy()
pvp.show()
