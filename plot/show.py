from argparse import ArgumentParser, BooleanOptionalAction
import pyvista as pv

parser = ArgumentParser()
parser.add_argument("fname", help="file that should be shown")
parser.add_argument("scalar", nargs="?", help="scalar to show", default=None)
parser.add_argument("-x", "--explode", help="show exploded mesh", action=BooleanOptionalAction)
parser.add_argument("-e", "--edges", help="show mesh edges", action=BooleanOptionalAction)
parser.add_argument("-c", "--cmap", help="colormap", default="coolwarm")
args = parser.parse_args()

mesh = pv.read(args.fname)
mesh = mesh.explode(1) if args.explode else mesh

pvp = pv.Plotter(lighting="three lights")
pvp.add_mesh(
    mesh,
    scalars=args.scalar,
    show_edges=args.edges,
    cmap=args.cmap,
    scalar_bar_args=dict(vertical=True),
)
pvp.show_bounds()
pvp.view_xy()
pvp.show()
