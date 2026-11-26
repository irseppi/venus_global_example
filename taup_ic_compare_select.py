"""
Record section of a few chosen distances (10, 60 and 175 deg), evenly
spaced on the y axis: solid inner core (V5-Th) vs. liquid core (V1-Tc),
with TauP travel times. The same plot as taup_ic_compare.py: each
station is scaled to the larger peak of its two runs, and that peak (m/s)
is written at the right of the trace.

Uses the run folders, period, phases and TauP models of
taup_ic_compare.py (change PERIOD there). Missing seismograms are
skipped.

Run with a Python that has obspy, numpy and matplotlib:
    conda activate venus
    python taup_ic_compare_select.py
"""
import numpy as np
from matplotlib import pyplot as plt

from taup_ic_compare import (EVENT_DEPTH_KM, NAMES, PERIOD, RUNS,
                             build_taup, read_ground_traces)

DISTANCES = [10, 60, 175]        # degrees
# P and S, then every single-pass phase that reaches the inner core:
# i = reflection off the ICB, I = P leg and J = S leg in the inner core;
# those that do not exist at a distance are left out
PHASES = ["P", "S",
          "PKiKP", "PKIKP", "PKIIKP", "PKJKP",
          "SKiKS", "SKIKS", "SKIIKS", "SKJKS",
          "PKiKS", "PKIKS", "PKJKS",
          "SKiKP", "SKIKP", "SKJKP"]
FILL = 0.9                       # each station's peak, in trace spacings
TICK_HALF = 0.4                  # half height of the TauP ticks, same unit
FIGURE = f"ic_compare_select_{PERIOD}.png"


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
        found = read_ground_traces(case_dir)
        traces[key] = {d: found[d] for d in DISTANCES if d in found}
        print(f"{case_dir}: {len(traces[key])} of {len(DISTANCES)} "
              "distances found")
        if traces[key]:
            models[key] = build_taup(case_dir, NAMES[key])
    dists = sorted(set(traces["solid"]) | set(traces["liquid"]))
    if not dists:
        print("no seismograms found")
        return

    colors = {ph: plt.cm.tab20(i % 20) for i, ph in enumerate(PHASES)}
    style = {"solid": dict(color="r", ls="-"),
             "liquid": dict(color="k", ls="--")}
    fig, ax = plt.subplots(figsize=(11, 5))
    labeled = set()
    t_max = 0.0
    # traces evenly spaced (one unit apart), whatever their distance
    for pos, dist in enumerate(dists):
        # one scale per station, shared by both runs: the larger peak of
        # the two fills FILL; that peak is written on the trace (in m/s)
        peak = max(np.abs(traces[key][dist][1]).max() for key in RUNS
                   if dist in traces[key])
        scale = FILL / max(peak, 1e-30)
        ax.text(1.005, pos, f"{peak:.1e}", transform=ax.get_yaxis_transform(),
                fontsize=7, va="center", ha="left", clip_on=False)
        for key in RUNS:
            if dist in traces[key]:
                t, y = traces[key][dist]
                ax.plot(t, y * scale + pos, lw=0.7, **style[key])
                t_max = max(t_max, t[-1])
        for key, ls in (("solid", "-"), ("liquid", ":")):
            if key not in models:
                continue
            for ph, ta in first_arrivals(models[key], dist).items():
                ax.vlines(ta, pos - TICK_HALF, pos + TICK_HALF,
                          colors=colors[ph], linestyles=ls, lw=0.8,
                          zorder=100, label=None if ph in labeled else ph)
                labeled.add(ph)

    ax.set_title("Vertical Velocity: Solid Inner Core (red) vs "
                 "Liquid Inner Core (black)\nTauP: Solid —, Liquid ···",
                 fontsize=10)
    ax.set_xlabel("Time (s)")
    ax.set_ylabel("Distance (degrees)")
    ax.set_yticks(range(len(dists)), [f"{d:g}" for d in dists])
    ax.set_xlim(0, t_max)
    # room for the largest swing above the top and below the bottom trace
    ax.set_ylim(-FILL - 0.1, len(dists) - 1 + FILL + 0.1)
    ax.text(1.005, 1.0, "peak\n(m/s)", transform=ax.transAxes, fontsize=7,
            va="bottom", ha="left")
    # legend in the order of PHASES
    handles, names = ax.get_legend_handles_labels()
    order = sorted(range(len(names)), key=lambda i: PHASES.index(names[i]))
    ax.legend([handles[i] for i in order], [names[i] for i in order],
              loc="upper left", fontsize=7, framealpha=0.9, ncol=2)
    fig.tight_layout()
    fig.savefig(FIGURE, dpi=150)
    print(f"wrote {FIGURE}")
    plt.show()


if __name__ == "__main__":
    main()
