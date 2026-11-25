#!/usr/bin/env python3
"""
Plot seismograms from the 'air' and 'noair' runs on top of each other.

Usage::

    python compare_runs.py runs/venus_T20 G030 G060 A030H050

Each station gets one figure with the two runs and their difference.
Needs numpy and matplotlib.
"""
import argparse
import glob
import os

import numpy as np
import matplotlib.pyplot as plt


def read_trace(run_dir, station, comp):
    """
    Read one SPECFEM2D ASCII seismogram.

    :type run_dir: str
    :param run_dir: case folder containing OUTPUT_FILES/
    :type station: str
    :param station: station name, e.g. G030
    :type comp: str
    :param comp: file suffix pattern, e.g. 'BXZ.semv' or 'PRE.semp'
    :rtype: tuple
    :return: (time, amplitude) arrays, or None if the file is missing
    """
    files = glob.glob(os.path.join(run_dir, "OUTPUT_FILES",
                                   f"*.{station}.{comp}"))
    if not files:
        return None
    data = np.loadtxt(files[0])
    return data[:, 0], data[:, 1]


def plot_station(exp_dir, station, cases, comp, out_dir, nyquist_freq=None):
    """
    Plot one station for all cases plus the air-minus-noair difference.

    :type exp_dir: str
    :param exp_dir: experiment folder holding one folder per case
    :type station: str
    :param station: station name
    :type cases: list
    :param cases: case names, first two are differenced
    :type comp: str
    :param comp: component suffix
    :type out_dir: str
    :param out_dir: where to save the PNG
    :type nyquist_freq: float or None
    :param nyquist_freq: optional low-pass cutoff in Hz, capped at each trace's Nyquist frequency
    """
    traces = {c: read_trace(os.path.join(exp_dir, c), station, comp)
              for c in cases}
    traces = {c: t for c, t in traces.items() if t is not None}
    if not traces:
        print(f"no {comp} data for {station}")
        return
    if nyquist_freq is not None:
        filtered_traces = {}
        for case, (t, a) in traces.items():
            if len(t) < 2:
                filtered_traces[case] = (t, a)
                continue
            dt = float(np.median(np.diff(t)))
            if dt <= 0:
                raise ValueError(f"non-increasing time samples for {case} at {station}")
            freqs = np.fft.rfftfreq(len(a), d=dt)
            cutoff = min(nyquist_freq, 0.5 / dt)
            spectrum = np.fft.rfft(a)
            spectrum[freqs > cutoff] = 0
            filtered_traces[case] = (t, np.fft.irfft(spectrum, n=len(a)))
        traces = filtered_traces
    fig, axes = plt.subplots(2, 1, figsize=(9, 5), sharex=True)
    for case, (t, a) in traces.items():
        axes[0].plot(t, a, lw=0.8, label=case)
    axes[0].set_ylabel(comp)
    axes[0].legend(loc="upper right")
    axes[0].set_title(station)
    if len(traces) >= 2:
        (t0, a0), (_, a1) = list(traces.values())[:2]
        axes[1].plot(t0, a0 - a1, lw=0.8, color="k")
        axes[1].set_ylabel(f"{cases[0]} - {cases[1]}")
    axes[1].set_xlabel("time (s)")
    fig.tight_layout()
    out = os.path.join(out_dir, f"{station}_{comp.replace('.', '_')}.png")
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"wrote {out}")


def main():
    """Command-line entry point."""
    ap = argparse.ArgumentParser(description="Compare air/noair runs")
    ap.add_argument("exp_dir", help="experiment folder, e.g. runs/venus_T20")
    ap.add_argument("stations", nargs="+",
                    help="stations to plot; omit to plot all stations")
    ap.add_argument("--cases", nargs="+", default=["air", "noair"])
    ap.add_argument("--nyquist-freq", nargs="?", const=float("inf"), type=float, default=0.05148756602499681*2,
                    help="optionally low-pass at this frequency in Hz; without a value, use each trace's Nyquist frequency")

    ap.add_argument("--comps", nargs="+",
                    default=["BXZ.semv", "BXX.semv", "PRE.semp", "HXZ.semv", "HXX.semv",])
    args = ap.parse_args()
    out_dir = os.path.join(args.exp_dir, "figures")
    os.makedirs(out_dir, exist_ok=True)
    for sta in args.stations:
        for comp in args.comps:
            print(args.nyquist_freq)
            plot_station(args.exp_dir, sta, args.cases, comp, out_dir,
                         args.nyquist_freq)


if __name__ == "__main__":
    main()

'''
python compare_runs.py runs/venus_T20 G010 G030 G060 A010H010 A010H050 \
  --nyquist-freq 0.05148756602499681

'''
