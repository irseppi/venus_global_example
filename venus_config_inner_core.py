"""
Settings for the solid vs. liquid inner core experiment (Venus model
V5-Th with a 2493 km inner core).

Run:  python make_venus_run.py venus_config_inner_core.py
This makes four run folders on the same mesh, every pair of air case
(CASES) and inner-core case (CORE_CASES):
    air_icsolid     Venus air, solid inner core from the V5-Th profile
    air_icliquid    Venus air, liquid control (see INNER_CORE_LIQUID)
    noair_icsolid   air density x 1e-3, solid inner core
    noair_icliquid  air density x 1e-3, liquid control
Run only the folders you need.
Nothing is run or submitted. In each run folder:
    run_venus_cluster.sh  ->  sbatch run_venus_cluster.sh
    run_venus_local.sh    ->  ./run_venus_local.sh > run.log 2>&1
"""

# --- SPECFEM2D builds ---------------------------------------------------------
CLUSTER_SPECFEM_BIN = "/home/irseppi/REPOS/specfem2d/bin"
LOCAL_SPECFEM_BIN = "/Users/irseppi/REPOSITORIES/specfem2d/bin"
CLUSTER_FC = "gfortran"
LOCAL_FC = "gfortran"

# --- output -------------------------------------------------------------------
OUTPUT_DIR = "runs"
EXPERIMENT = "venus_IC2500_T60"

# --- interior model -----------------------------------------------------------
# Made from data/models/V5-Th_ICB_2500.out with
#   python make_venus_run.py --convert-profile data/models/V5-Th_ICB_2500.out
# CMB 3183.3 km, inner-core boundary 2493.4 km (vs = 3.5 km/s inside).
INTERIOR_TABLE = "data/venus_interior_V5-Th_ICB_2500.txt"
INNER_CORE_RADIUS_KM = "table"   # take the inner core from the table

# case name -> True for a solid inner core, False for the liquid control.
# Every combination of CASES and CORE_CASES becomes a run folder.
CORE_CASES = {"icsolid": True, "icliquid": False}
# Liquid control: "outer_core" = outer core at the boundary held constant
# to the center (no inner core at all); "same_bulk" = same density and
# bulk modulus as the solid inner core, shear modulus 0 (only rigidity
# changes)
INNER_CORE_LIQUID = "outer_core"

# --- resolution ---------------------------------------------------------------
MIN_PERIOD_S = 60.0
POINTS_PER_WAVELENGTH = 5.0
MULTIPLICATION_FACTOR = None
INTERIOR_ASPECT_RATIO = 1.5    # interior elements ~1.5x taller than wide

# --- atmosphere ---------------------------------------------------------------
# "air" = Venus atmosphere; "noair" = air density x 1e-3 (same mesh).
# The sponge needs a SPECFEM2D build with
# fortran/axisym_viscoacoustic_fix.patch (the laptop build has it).
AIR_TOP_KM = 100.0
SPONGE = True
SPONGE_BOTTOM_KM = 70.0
SPONGE_Q_MIN = 3.0
CASES = {"air": 1.0, "noair": 1.0e-3}

# --- source -------------------------------------------------------------------
SOURCE_TYPE = "force"
SOURCE_DEPTH_KM = 10.0
SIM_DURATION_S = 3000.0
DT_OVERRIDE = None

# --- stations (angles in degrees from the source) -----------------------------
# Inner-core phases (PKIKP, PKiKP) need the far side of the planet
SEISMOTYPE = "2,4"
GROUND_STATION_DEG = list(range(5, 180, 5))     # 5, 10, ..., 175
AIR_STATION_DEG = list(range(5, 180, 5))        # above each ground station
AIR_STATION_ALT_KM = [1, 10, 30, 50, 60]        # keep below SPONGE_BOTTOM_KM

# --- cluster ------------------------------------------------------------------
NPROC = 168
SLURM = dict(
    job_name="venus_ic",
    nodes=5,
    ntasks=168,
    partition="t2standard",
    time="12:00:00",
    output="log.txt",
)

# --- laptop / workstation -----------------------------------------------------
LOCAL_NPROC = 4
LOCAL_MPIRUN = "mpirun"
