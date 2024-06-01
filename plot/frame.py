import os
import tqdm
import multiprocessing
import pyvista as pv
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from pypdf import PdfWriter
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument("fnames", nargs="+", help="vtkhdf files that should be plotted")
parser.add_argument("-s", "--scalar", help="the scalar to plot", default="density")
parser.add_argument("--xmin", help="minimum x value for plot", type=float, default=-np.inf)
parser.add_argument("--xmax", help="maximum x value for plot", type=float, default=np.inf)
parser.add_argument("--ymin", help="minimum y value for plot", type=float, default=-np.inf)
parser.add_argument("--ymax", help="maximum y value for plot", type=float, default=np.inf)
parser.add_argument("--cmap", help="the colormap to use", default="coolwarm")
parser.add_argument("--vmin", help="minimum value of scalar", type=float)
parser.add_argument("--vmax", help="maximum value of scalar", type=float)
args = parser.parse_args()


def plot(fname):
    fig, ax = plt.subplots(num=1, clear=True)

    pvp = pv.Plotter(off_screen=True, lighting="none")
    mesh = pv.read(fname, force_ext=".hdf")

    if len(mesh[args.scalar].shape) > 1:
        vmin = np.linalg.norm(mesh[args.scalar], axis=1).min()
        vmax = np.linalg.norm(mesh[args.scalar], axis=1).max()
    else:
        vmin = mesh[args.scalar].min()
        vmax = mesh[args.scalar].max()

    if args.vmin:
        vmin = args.vmin
    if args.vmax:
        vmax = args.vmax
    clim = [vmin, vmax]

    bounds = [args.xmin, args.xmax, args.ymin, args.ymax, -np.inf, np.inf]
    mesh = mesh.clip_box(bounds, invert=False)
    pvp.add_mesh(mesh, scalars=args.scalar, clim=clim, cmap=args.cmap, show_scalar_bar=False)

    pvp.camera.tight()
    img = pvp.screenshot(None, transparent_background=True, return_img=True, scale=2)

    extent = mesh.bounds[:4]

    im = ax.imshow(img, cmap=args.cmap, extent=extent, vmin=vmin, vmax=vmax)

    div = make_axes_locatable(ax)
    cax = div.append_axes("right", size=0.15, pad=0.1)
    cb = plt.colorbar(im, cax=cax)

    ax.set_xlabel("$x$")
    ax.set_ylabel("$y$")
    ax.set_title(f"time = {mesh["TimeValue"][0]:g}")
    cb.set_label(args.scalar)
    ax.set_xlim(xmin=max(args.xmin, extent[0]), xmax=min(args.xmax, extent[1]))
    ax.set_ylim(ymin=max(args.ymin, extent[2]), ymax=min(args.ymax, extent[3]))

    fig.tight_layout()
    fig.savefig(fname.replace(".vtkhdf", ".pdf"), dpi=300, bbox_inches="tight")


with multiprocessing.Pool(6) as pool:
    _ = list(tqdm.tqdm(pool.imap(plot, args.fnames), total=len(args.fnames)))

pdfs = sorted(list(map(lambda x: x.replace(".vtkhdf", ".pdf"), args.fnames)))
merger = PdfWriter()
for pdf in pdfs:
    merger.append(pdf)
merger.write(f"{pdfs[0].split('_')[0]}_{args.scalar.replace(" ", "_")}.pdf")
merger.close()
for pdf in pdfs:
    os.remove(pdf)
