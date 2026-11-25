"""
Plot the ground BXX velocity, the pressure 60 km up, and their difference
at the station farthest from the source (175 deg, the last station before
the antipode), for the solid inner-core run with atmosphere.

BXX is SPECFEM2D's x component (Cartesian, perpendicular to the symmetry
axis). Near the antipode that is almost the local horizontal: at 175 deg
it is 99.6% horizontal and 8.7% vertical motion.

Third panel, ground minus air:
    default               each trace divided by its own peak (-1 to 1),
                          so it compares waveform shape and timing
    --pressure-as-velocity  pressure turned into air particle velocity,
                          v = p / (rho c), with rho and c at ALT_KM from
                          the run's model table; difference in m/s. The
                          air velocity this gives is along the wave's
                          path (mostly vertical), while BXX is horizontal.

Vertical lines are TauP arrival times of the inner-core phases, from a
model built from the run's own DATA/venus_model.txt (see
taup_ic_compare.py). They are ground arrival times on every panel; the
pressure at altitude follows later by the travel time up through the air.

Run with a Python that has obspy, numpy and matplotlib:
    conda activate venus
    python plot_far_ground_vs_air.py
    python plot_far_ground_vs_air.py --pressure-as-velocity
"""
import argparse
import os

import numpy as np
from matplotlib import pyplot as plt

from taup_ic_compare import build_taup

RUN_DIR = "runs/venus_IC2500_T60/air_icsolid"
DIST_DEG = 10
ALT_KM = 60
EVENT_DEPTH_KM = 10.0
# inner-core phases; those that do not exist at DIST_DEG are left out
PHASES = ["PKiKP", "PKIKP", "PKIIKP", "SKiKS", "SKIKS", "PKJKP", "SKJKS"]
FIGURE = f"far_ground_vs_air_{DIST_DEG}deg.png"


def read_sem(path):
    """
    Read a SPECFEM2D ASCII seismogram.

    :type path: str
    :param path: seismogram file
    :rtype: tuple of numpy.ndarray
    :return: (time_s, data)
    """
    a = np.loadtxt(path)
    return a[:, 0], a[:, 1]


def air_impedance(case_dir, alt_km):
    """
    Air density and sound speed at an altitude, from the run's model
    table (region 4 rows of DATA/venus_model.txt).

    :type case_dir: str
    :param case_dir: run folder
    :type alt_km: float
    :param alt_km: altitude above the ground [km]
    :rtype: tuple of float
    :return: (rho [kg/m3], c [m/s])
    """
    with open(os.path.join(case_dir, "DATA", "venus_model.txt")) as f:
        rows = [ln.split() for ln in f
                if ln.strip() and not ln.startswith("#")][1:]
    solid = [float(p[1]) for p in rows if p[0] != "4"]
    air = sorted((float(p[1]), float(p[2]), float(p[3])) for p in rows
                 if p[0] == "4")
    r = max(solid) + alt_km * 1e3
    rr = [a[0] for a in air]
    return (float(np.interp(r, rr, [a[1] for a in air])),
            float(np.interp(r, rr, [a[2] for a in air])))


def main(argv=None):
    """
    Plot ground BXZ, air pressure and their difference.

    :type argv: list
    :param argv: command-line arguments (defaults to sys.argv[1:])
    """
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--pressure-as-velocity", action="store_true",
                    help="difference in m/s, with the pressure converted "
                         "to air particle velocity p/(rho c) (default: "
                         "both traces normalized to their peaks)")
    args = ap.parse_args(argv)

    out = os.path.join(RUN_DIR, "OUTPUT_FILES")
    t_g, ground = read_sem(os.path.join(
        out, f"VE.G{DIST_DEG:03d}.BXZ.semv"))
    t_a, press = read_sem(os.path.join(
        out, f"VE.A{DIST_DEG:03d}H{ALT_KM:03d}.PRE.semp"))
    air = np.interp(t_g, t_a, press)     # same time samples as the ground

    if args.pressure_as_velocity:
        rho, c = air_impedance(RUN_DIR, ALT_KM)
        diff = ground - air / (rho * c)
        diff_label = "BXZ $-$ p/($\\rho$c) (m/s)"
        diff_title = (f"ground $-$ air particle velocity, $\\rho$ = "
                      f"{rho:.3g} kg/m$^3$, c = {c:.0f} m/s at {ALT_KM} km")
        suffix = "_velocity"
    else:
        diff = ground / np.abs(ground).max() - air / np.abs(air).max()
        diff_label = "normalized difference"
        diff_title = "ground $-$ air, each divided by its own peak"
        suffix = "_normalized"

    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(9, 8.5), sharex=True)
    ax1.plot(t_g, ground, color="k", lw=0.8)
    ax1.set_ylabel("Ground velocity BXX (m/s)")
    ax1.set_title(f"Solid inner core (V5-Th), {DIST_DEG}$^\\circ$ from the "
                  "source")
    ax2.plot(t_a, press, color="b", lw=0.8)
    ax2.set_ylabel(f"Pressure at {ALT_KM} km (Pa)")
    ax3.plot(t_g, diff, color="m", lw=0.8)
    ax3.set_ylabel(diff_label)
    ax3.set_title(diff_title, fontsize=10)
    ax3.set_xlabel("Time (s)")
    ax3.set_xlim(0, max(t_g[-1], t_a[-1]))

    # earliest TauP arrival of each inner-core phase
    model = build_taup(RUN_DIR, "V5-Th_IC2500")
    arrivals = {}
    for a in model.get_travel_times(source_depth_in_km=EVENT_DEPTH_KM,
                                    distance_in_degree=DIST_DEG,
                                    phase_list=PHASES):
        arrivals[a.name] = min(a.time, arrivals.get(a.name, a.time))
    colors = {ph: plt.cm.tab10(i) for i, ph in enumerate(PHASES)}
    for ph, ta in sorted(arrivals.items(), key=lambda kv: kv[1]):
        for ax in (ax1, ax2, ax3):
            ax.axvline(ta, color=colors[ph], lw=1.2, ls="--",
                       label=f"{ph} ({ta:.0f} s)" if ax is ax1 else None)
    ax1.legend(loc="upper left", fontsize=8, framealpha=0.9,
               title="TauP inner-core phases", title_fontsize=8)
    for ax in (ax1, ax2, ax3):
        ax.grid(alpha=0.3)
    fig.tight_layout()
    name = FIGURE.replace(".png", suffix + ".png")
    fig.savefig(name, dpi=150)
    print(f"wrote {name}")
    plt.show()


if __name__ == "__main__":
    main()
