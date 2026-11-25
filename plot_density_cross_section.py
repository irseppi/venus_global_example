"""
Plot density cross sections (half disks, like a SPECFEM2D AXISYM mesh) of
the two Venus models as provided, including their atmospheres:

    data/models/V5-Th_ICB_2500.out   solid inner core
    data/models/V1_Steinberger_atm.out   V1-Tc (liquid core); same
                                         interior as V1-Tc.out plus the
                                         atmosphere, which V1-Tc.out lacks

The profiles are read straight from the .out files (first line a row
count, then radius_km T P rho_kg_m3 vp_km_s vs_km_s from the center up),
not from the SPECFEM run tables.

    python plot_density_cross_section.py
"""
import numpy as np
from matplotlib import pyplot as plt

MODELS = {"Solid Inner Core": "data/models/V5-Th_ICB_2500.out",
          "Liquid Inner Core": "data/models/V1_Steinberger_atm.out"}
CMAP = "turbo"
N_RADIUS = 3000          # radial cells (about 2 km each)
N_ANGLE = 720            # angular cells over the half disk
FIGURE = "density_cross_section.png"


def read_profile(path):
    """
    Read radius and density from a model .out file.

    :type path: str
    :param path: model file
    :rtype: tuple of numpy.ndarray
    :return: (radius_km, density_kg_m3), sorted by radius, with both
        values kept at each discontinuity
    """
    a = np.loadtxt(path, skiprows=1, usecols=(0, 3))
    order = np.argsort(a[:, 0], kind="stable")
    return a[order, 0], a[order, 1]


def main():
    """Draw both half-disk density cross sections with one colorbar."""
    profiles = {name: read_profile(p) for name, p in MODELS.items()}
    vmax = max(rho.max() for _, rho in profiles.values())
    fig, axes = plt.subplots(1, 2, figsize=(9, 7.5), facecolor="white")
    for ax, (name, (r, rho)) in zip(axes, profiles.items()):
        r_edges = np.linspace(0.0, r.max(), N_RADIUS + 1)
        r_mid = 0.5 * (r_edges[1:] + r_edges[:-1])
        rho_mid = np.interp(r_mid, r, rho)
        th = np.linspace(-np.pi / 2, np.pi / 2, N_ANGLE + 1)
        rr, tt = np.meshgrid(r_edges, th, indexing="ij")
        values = np.repeat(rho_mid[:, None], N_ANGLE, axis=1)
        mesh = ax.pcolormesh(rr * np.cos(tt), rr * np.sin(tt), values,
                             cmap=CMAP, vmin=0.0, vmax=vmax,
                             shading="flat", rasterized=True)
        ax.set_aspect("equal")
        ax.set_xlim(-0.03 * r.max(), 1.03 * r.max())
        ax.set_ylim(-1.03 * r.max(), 1.03 * r.max())
        ax.axis("off")
        ax.set_title(name, fontsize=11)
    cbar = fig.colorbar(mesh, ax=axes, shrink=0.7, pad=0.02)
    cbar.set_label("Density (kg/m$^3$)")
    fig.savefig(FIGURE, dpi=200, facecolor="white")
    print(f"wrote {FIGURE}")
    plt.show()


if __name__ == "__main__":
    main()
