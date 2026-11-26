"""
Plot the ground vertical velocity (local up, rotated from ?XX and ?XZ)
and the pressure ALT_KM up at one or more station distances, for the
solid inner-core run with atmosphere. Each distance gets its own figure.

Vertical lines are TauP arrival times of the inner-core phases, from a
model built from the run's own DATA/venus_model.txt (see
taup_ic_compare.py). TauP's model stops at the ground, so on the
pressure panel each arrival is moved later by the time its ray takes to
climb from the ground to ALT_KM through the run's air (region 4 of the
model table), using the arrival's own ray parameter. The rays are almost
vertical in the slow air, so this is close to the vertical travel time;
the small sideways drift of the ray in the air is ignored.

Run with a Python that has obspy, numpy and matplotlib:
    conda activate venus
    python plot_far_ground_vs_air.py              # distances in DISTANCES
    python plot_far_ground_vs_air.py 10 60 175    # these distances instead
"""
import argparse
import glob
import os

import numpy as np
from matplotlib import pyplot as plt

from taup_ic_compare import build_taup, read_local_vertical

PERIOD = "T20"                   # T20 or T60: picks the run folder
RUN_DIR = f"runs/venus_IC2500_{PERIOD}/air_icsolid"
DISTANCES = [10, 60, 175]        # degrees; one figure each
ALT_KM = 60
EVENT_DEPTH_KM = 10.0
# inner-core phases; those that do not exist at a distance are left out
PHASES = ["PKiKP", "PKIKP", "PKIIKP", "SKiKS", "SKIKS", "PKJKP", "SKJKS"]


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


def read_air(case_dir):
    """
    Air sound speed against radius, from the run's model table (region 4
    rows of DATA/venus_model.txt).

    :type case_dir: str
    :param case_dir: run folder
    :rtype: tuple
    :return: (ground radius [m], air radii [m], sound speeds [m/s])
    """
    with open(os.path.join(case_dir, "DATA", "venus_model.txt")) as f:
        rows = [ln.split() for ln in f
                if ln.strip() and not ln.startswith("#")][1:]
    r_ground = max(float(p[1]) for p in rows if p[0] != "4")
    air = sorted((float(p[1]), float(p[3])) for p in rows if p[0] == "4")
    return r_ground, np.array([a[0] for a in air]), np.array(
        [a[1] for a in air])


def air_delay(air, alt_km, ray_param):
    """
    Travel time of a ray from the ground up to an altitude through the
    air, t = integral of eta^2 / (r sqrt(eta^2 - p^2)) dr with eta = r/c.

    :type air: tuple
    :param air: output of read_air
    :type alt_km: float
    :param alt_km: altitude [km]
    :type ray_param: float
    :param ray_param: TauP ray parameter [s/rad]
    :rtype: float
    :return: delay [s]
    """
    r_ground, rr, cc = air
    r = np.linspace(r_ground, r_ground + alt_km * 1e3, 2001)
    eta = r / np.interp(r, rr, cc)
    return float(np.trapz(eta ** 2 / (r * np.sqrt(eta ** 2 - ray_param ** 2)),
                          r))


def plot_distance(model, air, dist_deg):
    """
    Plot the ground vertical velocity and the air pressure at one
    distance and save the figure.

    :type model: obspy.taup.TauPyModel
    :param model: TauP model of the run
    :type air: tuple
    :param air: air profile, output of read_air
    :type dist_deg: int
    :param dist_deg: station distance from the source [deg]
    :rtype: bool
    :return: False if the station's seismograms are missing
    """
    out = os.path.join(RUN_DIR, "OUTPUT_FILES")
    # the band letter (B, H, ...) depends on DT, so match any of them
    found = glob.glob(os.path.join(out, f"VE.G{dist_deg:03d}.?XZ.semv"))
    air_file = os.path.join(out, f"VE.A{dist_deg:03d}H{ALT_KM:03d}.PRE.semp")
    vertical = read_local_vertical(found[0], dist_deg) if found else None
    if vertical is None or not os.path.exists(air_file):
        print(f"{dist_deg} deg: no ground or {ALT_KM} km seismogram in "
              f"{out}, skipped")
        return False
    t_g, ground = vertical
    t_a, press = read_sem(air_file)

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(9, 6), sharex=True)
    ax1.plot(t_g, ground, color="k", lw=0.8)
    ax1.set_ylabel("Vertical Ground Velocity (m/s)")
    ax1.set_title(f"Solid Inner Core Model, {dist_deg}$^\\circ$ "
                  "from the source")
    ax2.plot(t_a, press, color="b", lw=0.8)
    ax2.set_ylabel(f"Pressure at {ALT_KM} km Altitude (Pa)")
    ax2.set_xlabel("Time (s)")
    ax2.set_xlim(0, max(t_g[-1], t_a[-1]))

    # earliest TauP arrival of each inner-core phase
    # (ground, air) times
    arrivals = {}
    for a in model.get_travel_times(source_depth_in_km=EVENT_DEPTH_KM,
                                    distance_in_degree=dist_deg,
                                    phase_list=PHASES):
        if a.name not in arrivals or a.time < arrivals[a.name][0]:
            arrivals[a.name] = (a.time,
                                a.time + air_delay(air, ALT_KM, a.ray_param))
    colors = {ph: plt.cm.tab10(i) for i, ph in enumerate(PHASES)}
    for ph, (tg, ta) in sorted(arrivals.items(), key=lambda kv: kv[1]):
        ax1.axvline(tg, color=colors[ph], lw=1.2, ls="--",
                    label=f"{ph} ({tg:.0f} s)")
        ax2.axvline(ta, color=colors[ph], lw=1.2, ls="--",
                    label=f"{ph} ({ta:.0f} s)")
    # near the source the waves arrive early (left), so put the legend right
    side = "upper right" if dist_deg < 40 else "upper left"
    ax1.legend(loc=side, fontsize=8, framealpha=0.9,
               title="TauP inner-core phases", title_fontsize=8)
    if arrivals:
        ax2.legend(loc=side, fontsize=8, framealpha=0.9,
                   title=f"at {ALT_KM} km", title_fontsize=8)
    for ax in (ax1, ax2):
        ax.grid(alpha=0.3)
    fig.tight_layout()
    name = f"far_ground_vs_air_{PERIOD}_{dist_deg}deg.png"
    fig.savefig(name, dpi=150)
    print(f"wrote {name}")
    return True


def main(argv=None):
    """
    Plot every requested distance, then show the figures.

    :type argv: list
    :param argv: command-line arguments (defaults to sys.argv[1:])
    """
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("distances", nargs="*", type=int, default=DISTANCES,
                    help=f"station distances in degrees (default: "
                         f"{DISTANCES})")
    args = ap.parse_args(argv)
    model = build_taup(RUN_DIR, "V5-Th_IC2500")
    air = read_air(RUN_DIR)
    for dist in args.distances:
        plot_distance(model, air, dist)
    plt.show()


if __name__ == "__main__":
    main()
