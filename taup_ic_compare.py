"""
Plot the T20 runs with atmosphere at the ground stations: V5-Th with a
solid inner core vs. V1-Tc with a liquid core, with TauP travel times for
both models.

V5-Th run in black, V1-Tc run in red (dashed), vertical velocity (?XZ),
each station scaled to the larger peak of the two runs there; that peak
(m/s) is written at the right of each trace. TauP ticks: V5-Th as solid
lines,
V1-Tc as dotted lines. The two models also differ in the mantle and CMB
radius (3183.3 vs 3352.75 km), not only in the inner core. The TauP
models are built from each run's own DATA/venus_model.txt, so they match
what SPECFEM used. Missing seismograms are skipped.

Run with a Python that has obspy, numpy and matplotlib:
    conda activate venus
    python taup_ic_compare.py
"""
import glob
import os

import numpy as np
from matplotlib import pyplot as plt
from obspy.taup import TauPyModel, taup_create

# run folders: V5-Th with a solid inner core, V1-Tc with a liquid core
PERIOD = "T20"                   # T20 or T60: picks the run folders
RUNS = {"solid": f"runs/venus_IC2500_{PERIOD}/air_icsolid",
        "liquid": f"runs/venus_V1Tc_{PERIOD}/air_V1Tc"}
NAMES = {"solid": "V5-Th_IC2500", "liquid": "V1-Tc"}
EVENT_DEPTH_KM = 10.0
PHASES = ["P", "PKiKP", "PKIKP", "S", "SKIKS"]
SPACING = 5.0                    # degrees between traces on the y axis
TAUP_DIR = "taup_models"
FIGURE = f"ic_compare_{PERIOD}.png"


def build_taup(case_dir, name):
    """
    Write a .tvel file (depth km, vp, vs km/s, rho g/cm3) from a run's
    DATA/venus_model.txt, without the air (region 4), and build a TauP
    model from it.

    :type case_dir: str
    :param case_dir: run folder
    :type name: str
    :param name: model name, also the file name
    :rtype: obspy.taup.TauPyModel
    :return: TauP model
    """
    with open(os.path.join(case_dir, "DATA", "venus_model.txt")) as f:
        lines = [ln.split() for ln in f
                 if ln.strip() and not ln.startswith("#")][1:]
    # (region, radius, rho, vp, vs); regions 1-3 get deeper with the flag
    rows = [(int(p[0]),) + tuple(float(v) for v in p[1:5]) for p in lines
            if p[0] != "4"]
    rows.sort(key=lambda r: (-r[1], r[0]))
    if rows[-1][1] > 0:                       # extend to the center
        rows.append(rows[-1][:1] + (0.0,) + rows[-1][2:])
    r_surf = rows[0][1]
    os.makedirs(TAUP_DIR, exist_ok=True)
    tvel = os.path.join(TAUP_DIR, name + ".tvel")
    with open(tvel, "w") as f:
        f.write(f"{name} - P\n{name} - S\n")
        last = None
        for _, r, rho, vp, vs in rows:
            line = (f"{(r_surf - r) / 1e3:10.4f} {vp / 1e3:10.4f} "
                    f"{vs / 1e3:10.4f} {rho / 1e3:10.4f}\n")
            if line != last:      # skip repeated rows (region boundaries
                f.write(line)     # with no jump, e.g. mid-mantle)
            last = line
    taup_create.build_taup_model(tvel, output_folder=TAUP_DIR)
    return TauPyModel(model=os.path.join(TAUP_DIR, name + ".npz"))


def read_ground_traces(case_dir):
    """
    Read the vertical velocity at the ground stations (VE.G<deg>).

    :type case_dir: str
    :param case_dir: run folder
    :rtype: dict
    :return: {distance_deg: (time, data)}; empty if nothing is there
    """
    out = {}
    pattern = os.path.join(case_dir, "OUTPUT_FILES", "VE.G*.?XZ.semv")
    for path in glob.glob(pattern):
        dist = float(os.path.basename(path).split(".")[1][1:])
        try:
            a = np.loadtxt(path, ndmin=2)
        except ValueError:
            continue
        if len(a) > 1:
            out[dist] = (a[:, 0], a[:, 1])
    return out


def first_arrivals(model, dist):
    """
    Earliest arrival time of each phase in PHASES.

    :type model: obspy.taup.TauPyModel
    :param model: TauP model
    :type dist: float
    :param dist: distance [deg]
    :rtype: dict
    :return: {phase name: time [s]}
    """
    out = {}
    for a in model.get_travel_times(source_depth_in_km=EVENT_DEPTH_KM,
                                    distance_in_degree=dist,
                                    phase_list=PHASES):
        out[a.name] = min(a.time, out.get(a.name, a.time))
    return out


def main():
    """Read both runs, build both TauP models and plot."""
    traces, models = {}, {}
    for key, case_dir in RUNS.items():
        traces[key] = read_ground_traces(case_dir)
        print(f"{case_dir}: {len(traces[key])} ground seismograms")
        if traces[key]:
            models[key] = build_taup(case_dir, NAMES[key])
    dists = sorted(set(traces["solid"]) | set(traces["liquid"]))
    if not dists:
        print("no seismograms found")
        return

    colors = {"P": "#1f77b4", "PKiKP": "#ff7f0e", "PKIKP": "#2ca02c",
              "S": "#9467bd", "SKIKS": "#17becf"}
    style = {"solid": dict(color="r", ls="-"),
             "liquid": dict(color="k", ls="--")}
    half = 0.4 * SPACING
    fig, ax = plt.subplots(figsize=(6, 8))
    labeled = set()
    t_max = 0.0
    for dist in dists:
        # one scale per station, shared by both runs: the larger peak of
        # the two fills 90% of the spacing; that peak is written on the
        # trace (in m/s) so amplitudes can be read off
        peak = max(np.abs(traces[key][dist][1]).max() for key in RUNS
                   if dist in traces[key])
        scale = 0.9 * SPACING / max(peak, 1e-30)
        ax.text(1.01, dist, f"{peak:.1e}", transform=ax.get_yaxis_transform(),
                fontsize=6, va="center", ha="left", clip_on=False)
        for key in RUNS:
            if dist in traces[key]:
                t, y = traces[key][dist]
                ax.plot(t, y * scale + dist, lw=0.7, **style[key])
                t_max = max(t_max, t[-1])
        for key, ls in (("solid", "-"), ("liquid", ":")):
            if key not in models:
                continue
            for ph, ta in first_arrivals(models[key], dist).items():
                ax.vlines(ta, dist - half, dist + half, colors=colors[ph],
                          linestyles=ls, lw=0.8, zorder=100,
                          label=None if ph in labeled else ph)
                labeled.add(ph)

    ax.set_title(f"Vertical Velocity: Solid Inner Core (red) vs "
                 "Liquid Inner Core (black)\nTauP: Solid —, Liquid ···",
                 fontsize=10)
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Distance (degrees)")
    ax.text(1.01, 1.0, "peak\n(m/s)", transform=ax.transAxes, fontsize=6,
            va="bottom", ha="left")
    ax.set_xlim(0, t_max)
    ax.set_ylim(0, 180)
    ax.legend(loc="lower right", fontsize=8, framealpha=0.9)
    fig.tight_layout()
    fig.savefig(FIGURE, dpi=150)
    print(f"wrote {FIGURE}")
    plt.show()


if __name__ == "__main__":
    main()
