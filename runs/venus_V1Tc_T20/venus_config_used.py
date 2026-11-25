"""
Settings for the 20 s V1-Tc liquid-core run: the companion of the 20 s
V5-Th solid inner-core run (venus_config_inner_core_T20.py).

V1-Tc has a fully liquid core and its CMB at 3352.75 km (V5-Th: 3183.3
km), so it needs its own mesh. Everything else (period, air, sponge,
source, stations, duration) matches venus_config_inner_core_T20.py, so the
seismograms can be compared station by station.

Run:  python make_venus_run.py venus_config_v1tc_T20.py
This makes one run folder:
    runs/venus_V1Tc_T20/air_V1Tc   Venus air, V1-Tc interior
Nothing is run or submitted. In the run folder:
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
EXPERIMENT = "venus_V1Tc_T20"

# --- interior model -----------------------------------------------------------
# Made from data/models/V1-Tc.out with
#   python make_venus_run.py --convert-profile data/models/V1-Tc.out
# CMB 3352.75 km, liquid core to the center.
INTERIOR_TABLE = "data/venus_interior_V1-Tc.txt"
# The mesh uses the same zone layout as the V5-Th run, with its inner-core
# zone at 2493.38 km. Here that radius is only a mesh boundary inside the
# liquid core: the material is V1-Tc's liquid core on both sides.
INNER_CORE_RADIUS_KM = 2493.3828
CORE_CASES = {"V1Tc": False}

# --- resolution ---------------------------------------------------------------
MIN_PERIOD_S = 20.0
POINTS_PER_WAVELENGTH = 5.0
MULTIPLICATION_FACTOR = None
INTERIOR_ASPECT_RATIO = 1.5

# --- atmosphere ---------------------------------------------------------------
AIR_TOP_KM = 100.0
SPONGE = True
SPONGE_BOTTOM_KM = 70.0
SPONGE_Q_MIN = 3.0
CASES = {"air": 1.0}

# --- source -------------------------------------------------------------------
SOURCE_TYPE = "force"
SOURCE_DEPTH_KM = 10.0
SIM_DURATION_S = 3000.0
DT_OVERRIDE = None

# --- stations (angles in degrees from the source) -----------------------------
SEISMOTYPE = "2,4"
GROUND_STATION_DEG = list(range(5, 180, 5))     # 5, 10, ..., 175
AIR_STATION_DEG = list(range(5, 180, 5))
AIR_STATION_ALT_KM = [1, 10, 30, 50, 60]

# --- cluster ------------------------------------------------------------------
NPROC = 168
SLURM = dict(
    job_name="venus_v1tc_T20",
    nodes=5,
    ntasks=168,
    partition="t2standard",
    time="12:00:00",
    output="log.txt",
)

# --- laptop / workstation -----------------------------------------------------
LOCAL_NPROC = 4
LOCAL_MPIRUN = "mpirun"
