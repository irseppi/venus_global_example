

import os
import glob
import numpy as np

from pysep import logger
from collections import defaultdict
from matplotlib import pyplot as plt
from obspy import UTCDateTime, Stream, Trace
from obspy.taup import TauPyModel, taup_create
from pysep.utils.cap_sac import append_sac_headers, append_sac_headers_cartesian
from obspy.taup import plot_travel_times
def read_sem(fid, origintime="1970-01-01T00:00:00", source=None, stations=None, 
             location="", precision=4, source_format="CMTSOLUTION"):
    """
    Specfem3D outputs seismograms to ASCII (.sem? or .sem.ascii) files.
    Converts SPECFEM synthetics into ObsPy Stream objects with the correct
    header information. If `source` and `stations` files are also provided,
    PySEP will write appropriate SAC headers to the underlying data.

    :type fid: str
    :param fid: path of the given ascii file
    :type origintime: obspy.UTCDateTime
    :param origintime: UTCDatetime object for the origintime of the event. If
        None given, defaults to dummy value of '1970-01-01T00:00:00'
    :type source: str
    :param source: optional SPECFEM source file (e.g., CMTSOLUTION, SOURCE)
        defining the event which generated the synthetics. Used to grab event
        information and append as SAC headers to the ObsPy Stream
    :type stations: str
    :param stations: optional STATIONS file defining the station locations for
        the SPECFEM generated synthetics, used to generate SAC headers
    :type location: str
    :param location: location value for a given station/component
    :type precision: int
    :param precision: dt precision determined by differencing two
        adjancent time steps in the underlying ascii text file.
    :rtype st: obspy.Stream.stream
    :return st: stream containing header and data info taken from ascii file
    """
    # This was tested up to SPECFEM3D Cartesian git version 6895e2f7
    try:
        times = np.loadtxt(fname=fid, usecols=0)
        data = np.loadtxt(fname=fid, usecols=1)

    # At some point in 2018, the Specfem developers changed how the ascii files
    # were formatted from two columns to comma separated values, and repeat
    # values represented as 2*value_float where value_float represents the data
    # value as a float
    except ValueError:
        times, data = [], []
        with open(fid, 'r') as f:
            lines = f.readlines()
        for line in lines:
            try:
                time_, data_ = line.strip().split(',')
            except ValueError:
                if "*" in line:
                    time_ = data_ = line.split('*')[-1]
                else:
                    raise ValueError
            times.append(float(time_))
            data.append(float(data_))

        times = np.array(times)
        data = np.array(data)

    # We assume that dt is constant after 'precision' decimal points
    delta = round(times[1] - times[0], precision)

    # Get metadata information from CMTSOLUTION and STATIONS files
    event = None
    if source is None:
        origintime = UTCDateTime(origintime)
    else:
        event = read_events_plus(source, format=source_format)[0]
        origintime = event.preferred_origin().time
        logger.info(f"reading origintime from event: {origintime}")

    # Honor that Specfem doesn't start exactly on 0 due to USER_T0
    origintime += times[0]

    # SPECFEM2D/SPECFEM3D_Cartesian style name format, e.g., NZ.BFZ.BXE.semd OR
    # SPECFEM3D_Globe style name format, e.g., TA.O20K.BXR.sem.ascii
    net, sta, cha, fmt, *_ = os.path.basename(fid).split(".")
    stats = {"network": net, "station": sta, "location": location,
             "channel": cha, "starttime": origintime, "npts": len(data),
             "delta": delta, "mseed": {"dataquality": 'D'}, "format": fmt
             }
    st = Stream([Trace(data=data, header=stats)])

    if event and stations:
        try:
            # `read_stations` will throw a ValueError for Cartesian coordinates
            inv = read_stations(stations)
            st = append_sac_headers(st, event, inv)
        except ValueError as e:
            # If Cartesian coordinate system, slightly different header approach
            st = append_sac_headers_cartesian(st, event, stations)
        # Broad catch here as this is an optional step that might not always
        # work or be possible
        except Exception as e:
            logger.warning(f"could not append SAC header to trace because {e}")

    return st

event_depth = 10 #km
cutoff = 0.05148756602499681 #from specfem

taup_create.build_taup_model("V1-Tc.tvel")
model = TauPyModel(model="V1-Tc")

component = 'vertical' # 'vertical' or 'horizontal' 
run_type = 60 # 20 or 60
phase_dict = defaultdict(list)
phase_colors = {}
sensor_type = 'both' # 'ground', 'air', or 'both'
if sensor_type == 'ground':
    trs_air = {}
    trs_noair = {}
elif sensor_type == 'air':
    trs_air = {}
elif sensor_type == 'both':
    trs_pressure = {}
    trs_ground = {}
if run_type == 60:
    folders = ['runs/venus_T60/noair/OUTPUT_FILES/', 
               'runs/venus_T60/air/OUTPUT_FILES/']
    if component == 'vertical':
        file_end = '.BXZ.semv'
    elif component == 'horizontal':
        file_end = '.BXX.semv'

elif run_type == 20:
    folders = ['runs/venus_T20/noair/OUTPUT_FILES/', 
               'runs/venus_T20/air/OUTPUT_FILES/']
    if component == 'vertical':
        file_end = '.HXZ.semv'
    elif component == 'horizontal':
        file_end = '.HXX.semv'

if sensor_type == 'ground':
    file_type = 'VE.G*'
elif sensor_type == 'air':
    folders = [folders[1]] # only air folder
    file_type = 'VE.A*'
    file_end = '30.PRE.semp'
elif sensor_type == 'both':
    folders = [folders[1]] # both folders


for folder in folders:
    fids = []
    if sensor_type == 'both':
        for file_end in ['.BXZ.semv', '30.PRE.semp']:
            if file_end == '.BXZ.semv':
                file_type = 'VE.G*'
            elif file_end == '30.PRE.semp':
                file_type = 'VE.A*'
            
            for fid in glob.glob(os.path.join(folder, file_type + file_end)):
                fids.append(fid)
    else:
        for fid in glob.glob(os.path.join(folder, file_type + file_end)):
            fids.append(fid)

    for jj, fid in enumerate(fids):
        print(fid)
        st = read_sem(fid, origintime="1970-01-01T00:00:00",
                    source=None, stations=None, location="", precision=4,
                    source_format="CMTSOLUTION")
                    
        tr = st[0]
        filename = os.path.basename(fid)

        if filename.startswith("VE.G"):
            print(filename.split('.')[1][1:])
            dist_degree = float(filename.split('.')[1][1:])
        elif filename.startswith("VE.A"):
            print(filename.split('.')[1][2:4])
            dist_degree = float(filename.split('.')[1][2:4])


        #tr.filter('lowpass', freq=cutoff, corners=4, zerophase=True)

        if 'noair' in folder:
            trs_noair[dist_degree] = tr.copy()
        elif 'air' in folder:
            if sensor_type == 'both':
                print(file_type)
                if filename.startswith("VE.G"):
                    print(dist_degree)
                    trs_ground[dist_degree] = tr.copy()
                elif filename.startswith("VE.A"):
                    trs_pressure[dist_degree] = tr.copy()
            else:
                trs_air[dist_degree] = tr.copy()

        arrival_s = model.get_travel_times(source_depth_in_km=event_depth,
                                            distance_in_degree=dist_degree,
                                            phase_list=['ttbasic'])#['P', 'S', 'PKP', 'SKS', 'PKiKP', 'PKIKP']) #'ttall']) 
        # Keep only the earliest (first) arrival for each phase at this distance.
        first_arrivals = {}
        for a in arrival_s:
            if a.name not in first_arrivals or a.time < first_arrivals[a.name].time:
                first_arrivals[a.name] = a
        plot_each = False
        if plot_each:
            plt.figure()
            plt.plot(tr.times(), tr.data, color='k', linewidth=0.7)
        for a in first_arrivals.values():
            phase_dict[a.name].append((dist_degree, a.time))
            phase_colors.setdefault(
                a.name, plt.cm.tab20(len(phase_colors) % 20)
            )
            if plot_each: 
                plt.axvline(a.time, linewidth=1.0,
                        color=phase_colors[a.name], label=a.name, zorder=100)
        if plot_each:
            plt.title(f"{filename} at {dist_degree} degrees")
            plt.xlabel("Time (s)")
            plt.ylabel("Amplitude")
            plt.legend(loc='upper right', fontsize=8, frameon=True, framealpha=0.9)
            plt.show()
            plt.close()
# determine normalization scale (so traces visually comparable)
spacing = 5 # vertical spacing between traces in degrees (in y units)
scale_by = 'individual' # 'individual', 'global' or set
if sensor_type != 'both' and scale_by == 'global':
    max_amp = max((np.nanmax(np.abs(trs_air[dist_degree].data)) for dist_degree in trs_air.keys()))
    amp_scale = 0.9 * spacing / max_amp  # scale so peaks occupy ~80% of spacing

colors = {phase: plt.cm.tab20(i % 20) for i, phase in enumerate(phase_dict.keys())}
if scale_by == 'set':
    amp_scale = 5000000000
fig, ax = plt.subplots(1, 1, figsize=(4, 5))
if sensor_type != 'both':
    for d, tr in trs_noair.items():
        t = tr.times()
        if scale_by == 'individual':
            amp_scale = 0.9 * spacing / max(np.abs(tr.data))  # scale so peaks occupy ~80% of spacing
        scaled = tr.data * amp_scale + d
        ax.plot(t, scaled, color='k', linewidth=0.8)

    for d, tr in trs_air.items():
        t = tr.times()
        if scale_by == 'individual':
            amp_scale = 0.9 * spacing / max(np.abs(tr.data))  # scale so peaks occupy ~80% of spacing
        scaled = tr.data * amp_scale + d
        ax.plot(t, scaled, color='b', linestyle='dotted', linewidth=0.7)

elif sensor_type == 'both':
    print(f"trs_ground: {trs_ground.keys()}")
    for d, tr in trs_ground.items():
        t = tr.times()
        if scale_by == 'individual':
            amp_scale = 0.9 * spacing / max(np.abs(tr.data))  # scale so peaks occupy ~80% of spacing
        scaled = tr.data * amp_scale + d
        ax.plot(t, scaled, color='k', linewidth=0.8)

    for d, tr in trs_pressure.items():
        t = tr.times()
        if scale_by == 'individual':
            amp_scale = 0.9 * spacing / max(np.abs(tr.data))  # scale so peaks occupy ~80% of spacing
        elif scale_by == 'set':
            amp_scale = 5000000000.000
        scaled = tr.data * amp_scale + d
        ax.plot(t, scaled, color='b', linestyle='dotted', linewidth=0.7)
arrival_half_height = 0.4 * 0.8 * spacing
for phase in phase_dict.keys():
    pts = phase_dict[phase]
    ys, xs = [d for d, t in pts], [t for d, t in pts]

    ax.vlines(xs, ymin=[y - arrival_half_height for y in ys],
            ymax=[y + arrival_half_height for y in ys],
            colors=colors[phase], linewidth=1.0, alpha=0.9,
            label=phase, zorder=100)

ax.set_ylabel(f"Distance (degrees)") #(x{amp_scale:.1e})")
ax.set_xlabel("Time (s)")
ax.legend(loc='upper right', fontsize=8, frameon=True, framealpha=0.9)
ax.set_xlim(0, max(tr.times()))
plt.show()
