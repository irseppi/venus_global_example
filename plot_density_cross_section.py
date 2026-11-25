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

A second figure is a flat (square) cross section of the atmosphere
(identical in both files) above the crust, with a blue log-scale density
colorbar. Distances are true to scale; the curvature is left out.

    python plot_density_cross_section.py
"""
import numpy as np
from matplotlib import pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.patches import Rectangle

MODELS = {"Solid Inner Core": "data/models/V5-Th_ICB_2500.out",
          "Liquid Inner Core": "data/models/V1_Steinberger_atm.out"}
CMAP = "turbo"
N_RADIUS = 3000          # radial cells (about 2 km each)
N_ANGLE = 720            # angular cells over the half disk
FIGURE = "density_cross_section.png"
ATM_MODEL = "data/models/V5-Th_ICB_2500.out"
ATM_CMAP = "Blues"
CRUST_KM = 20.0          # drawn crust depth (the 2950 kg/m3 layer
CRUST_RHO = 2950.0       # goes to 48-60 km in both models)
CRUST_COLOR = "#c8a27a"
ATM_FIGURE = "density_cross_section_atmosphere.png"


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


def read_atmosphere(path):
    """
    Read the atmosphere rows (vs = 0 above the solid surface) of a model
    .out file.

    :type path: str
    :param path: model file
    :rtype: tuple of numpy.ndarray
    :return: (altitude_km, density_kg_m3), sorted by altitude
    """
    a = np.loadtxt(path, skiprows=1, usecols=(0, 3, 5))
    r_ground = a[a[:, 2] > 0, 0].max()
    air = a[(a[:, 0] >= r_ground) & (a[:, 2] == 0)]
    order = np.argsort(air[:, 0], kind="stable")
    return air[order, 0] - r_ground, air[order, 1]


def plot_atmosphere():
    """
    Draw the atmosphere's density as a flat (square) cross section above
    a labeled crust; the profile is 1-D, so it is the same across.
    """
    alt, rho = read_atmosphere(ATM_MODEL)
    top = alt.max()
    a_edges = np.linspace(0.0, top, 401)
    a_mid = 0.5 * (a_edges[1:] + a_edges[:-1])
    # log interpolation: density falls by about six decades
    rho_mid = np.exp(np.interp(a_mid, alt, np.log(rho)))
    width = top + CRUST_KM                       # square: width = height
    x_edges = np.array([0.0, width])

    fig, ax = plt.subplots(figsize=(7, 6), facecolor="white")
    mesh = ax.pcolormesh(x_edges, a_edges, rho_mid[:, None], cmap=ATM_CMAP,
                         norm=LogNorm(rho_mid.min(), rho_mid.max()),
                         shading="flat", rasterized=True)
    ax.add_patch(Rectangle((0, -CRUST_KM), width, CRUST_KM, fc=CRUST_COLOR,
                           ec="none"))
    ax.axhline(0, color="k", lw=1.0)
    ax.text(width / 2, -CRUST_KM / 2,
            f"Crust ($\\rho$ = {CRUST_RHO:.0f} kg/m$^3$)",
            ha="center", va="center", fontsize=11)
    ax.set_aspect("equal")
    ax.set_xlim(0, width)
    ax.set_ylim(-CRUST_KM, top)
    ax.set_yticks(np.arange(-20 * (CRUST_KM // 20), top + 1, 20))
    ax.set_xlabel("Horizontal distance (km)")
    ax.set_ylabel("Altitude (km)")
    ax.set_title("Venus Atmosphere", fontsize=11)
    cbar = fig.colorbar(mesh, ax=ax, shrink=0.8, pad=0.03)
    cbar.set_label("Atmosphere density (kg/m$^3$)")
    fig.savefig(ATM_FIGURE, dpi=200, facecolor="white",
                bbox_inches="tight")
    print(f"wrote {ATM_FIGURE}")


def main():
    """Draw both half-disk density cross sections, then the atmosphere."""
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
    plot_atmosphere()
    plt.show()


if __name__ == "__main__":
    main()
