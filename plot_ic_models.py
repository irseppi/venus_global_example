"""
Plot P velocity, S velocity and density against radius for V5-Th (solid
inner core) and V1-Tc (liquid core), without the atmosphere.

The models are read from the T60 run folders' DATA/venus_model.txt, so
they are exactly what SPECFEM uses. V5-Th in black, V1-Tc in red
(dashed).

    python plot_ic_models.py
"""
import os

from matplotlib import pyplot as plt

RUNS = {"V5-Th, solid inner core": "runs/venus_IC2500_T60/air_icsolid",
        "V1-Tc, liquid core": "runs/venus_V1Tc_T60/air_V1Tc"}
FIGURE = "ic_models.png"


def read_model(case_dir):
    """
    Read a run's DATA/venus_model.txt without the air (region 4).

    :type case_dir: str
    :param case_dir: run folder
    :rtype: tuple of list
    :return: (radius_km, vp_km_s, vs_km_s, rho_g_cm3), ordered from the
        surface down, with both values at each discontinuity
    """
    with open(os.path.join(case_dir, "DATA", "venus_model.txt")) as f:
        lines = [ln.split() for ln in f
                 if ln.strip() and not ln.startswith("#")][1:]
    rows = [(int(p[0]),) + tuple(float(v) for v in p[1:5]) for p in lines
            if p[0] != "4"]
    # regions 1-3 (mantle, outer core, inner core) get deeper with the flag
    rows.sort(key=lambda r: (-r[1], r[0]))
    return ([r[1] / 1e3 for r in rows], [r[3] / 1e3 for r in rows],
            [r[4] / 1e3 for r in rows], [r[2] / 1e3 for r in rows])


def main():
    """Plot vp, vs and density of both models against radius."""
    styles = [dict(color="k", ls="-", lw=1.5),
              dict(color="r", ls="--", lw=1.5)]
    fig, axes = plt.subplots(1, 3, figsize=(10, 6), sharey=True)
    for (label, case_dir), style in zip(RUNS.items(), styles):
        r, vp, vs, rho = read_model(case_dir)
        axes[0].plot(vp, r, label=label, **style)
        axes[1].plot(vs, r, label=label, **style)
        axes[2].plot(rho, r, label=label, **style)
    for ax, name in zip(axes, ("P velocity (km/s)", "S velocity (km/s)",
                               "Density (g/cm$^3$)")):
        ax.set_xlabel(name)
        ax.grid(alpha=0.3)
    axes[0].set_ylabel("Radius (km)")
    axes[0].set_ylim(0, None)
    axes[0].legend(loc="lower left", fontsize=9)
    fig.suptitle("Venus interior models")
    fig.tight_layout()
    fig.savefig(FIGURE, dpi=150)
    print(f"wrote {FIGURE}")
    plt.show()


if __name__ == "__main__":
    main()
