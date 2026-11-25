# CLAUDE.md: Venus global axisymmetric SPECFEM2D example

Context for Claude sessions that continue this work. Read this first, then
`README.md` (user-facing steps). Last updated 23 Nov 2026.

## Project

NASA-funded study preparing for a landed seismometer on Venus: how does
Venus's thick atmosphere affect what a surface seismometer records? This
folder builds global **axisymmetric (AXISYM) SPECFEM2D** simulations of
Venus with a fluid atmosphere (0–100 km), a solid mantle, and a liquid
core. The core science run is a pair of simulations on the same mesh,
**with vs. without atmosphere**, compared at ground stations and at air
stations (1–60 km, balloon altitudes).

Near-term goal: results for an AGU26 poster (December 2026). Planned
next: solid vs. liquid inner core (generator support added 23 Nov 2026;
see "Inner-core layout" below).

## Coding conventions (from the PI's global CLAUDE.md)

- Maximum line width 80 characters.
- Python follows PEP 8.
- Every function and module docstring has a short description; function
  docstrings list parameters with types in reST style
  (`:type x: float` / `:param x: ...`, plus `:rtype:` / `:return:`).
- `make_venus_run.py` uses only the Python standard library (Python ≥
  3.7). Keep it that way so it runs on the cluster without installs.
  `compare_runs.py` may use numpy and matplotlib.

## Ground rules

- The generator **only writes files**. It never runs the mesher or the
  solver and never submits jobs; the user runs them.
- Do not modify the user's SPECFEM2D checkout or run anything on the
  cluster unless asked. Test with a scratch copy of SPECFEM2D.
- Changes to SPECFEM2D source ship as files in `fortran/` with README
  instructions, never as silent edits.

## Files

| File | Role |
|---|---|
| `venus_config.py` | The only file users normally edit (Python module, UPPERCASE settings) |
| `venus_config_inner_core.py` | Solid vs. liquid inner-core experiment (V5-Th model) |
| `data/models/V5-Th_ICB_2500.out` | V5-Th profile (r_km T P rho vp vs) with a solid inner core |
| `data/venus_interior_V5-Th_ICB_2500.txt` | It, converted with `--convert-profile` (mantle 1, core 3, Q 9999) |
| `make_venus_run.py` | Generator: config → one run folder per case |
| `compare_runs.py` | Overlays two runs' seismograms plus their difference |
| `fortran/define_external_model.f90` | Replacement SPECFEM2D model routine: reads `DATA/venus_model.txt` |
| `fortran/axisym_viscoacoustic_fix.patch` | Bug fix for fluid attenuation in AXISYM (needed for the sponge) |
| `data/venus_interior.txt` | Solid planet and core, region-tagged, SI units |
| `data/venus_air_profile.txt` | Approximate VIRA atmosphere (alt, T, P, rho, c) |
| `templates/create_mesh_Venus.F90` | Stand-alone Fortran mesher (parameters patched by the generator) |
| `templates/Par_file`, `templates/SOURCE` | From SPECFEM2D `EXAMPLES/applications/axisymmetric_examples/global_Earth_ak135f/DATA` |
| `README.md` | User-facing step-by-step guide |

## How a run is built

```text
venus_config.py
   │
make_venus_run.py ──► runs/<EXPERIMENT>/<case>/
   ├─ DATA/venus_model.txt      interior rows + air rows (+ sponge Q, density scale)
   ├─ DATA/Par_file, SOURCE, STATIONS   (key = value edits of the templates)
   ├─ MESH/Nummaterial_velocity_file_NO_MUD  (region → fluid/solid)
   ├─ create_mesh_venus.F90      (mesher with radii and element counts patched)
   ├─ run_venus_cluster.sh       (SLURM; sbatch from inside the folder)
   ├─ run_venus_local.sh         (laptop; run from inside the folder)
   └─ run_summary.txt            (resolution, DT, NSTEP, echo time, cost)
```

Each run script compiles and runs the mesher, rewrites `NPROC` in the
Par_file to its own process count, links `xmeshfem2D`/`xspecfem2D` from
`CLUSTER_SPECFEM_BIN` or `LOCAL_SPECFEM_BIN`, and stops early if the solver
was not built with the Venus model reader (it greps the binary for
"Venus 1-D table"). The scripts must not `cd "$(dirname "$0")"`, because
SLURM runs a spooled copy of the script.

### Mesh and region flags

The mesher is the SPECFEM2D global AK135F mesher with the Earth layers
shifted down one slot to make room for air. It has five fixed radial
zones; the flags it writes and their current meaning:

| Mesher zone (Earth name) | Flag | Current Venus meaning | Fluid? |
|---|---|---|---|
| d670 → surface (`R670` → `R_EARTH`) | 4 | Air, 0–100 km | yes |
| CMB → d670 (`RCMB` → `R670`) | 1 | Upper mantle, surface to 738 km depth | no |
| ICB → doubling → CMB (`RICB` → `RCMB`) | 2 | Lower mantle | no |
| Cube → ICB, plus central cube (< `RICB`) | 3 | Whole core | yes |

Mesher radii are written from the model table: `R670` = 6051.9 km
(surface), `RCMB` = 5313.568 km (mid-mantle boundary), `RICB` = 3352.753
km (core–mantle boundary), `R_EARTH` = surface + `AIR_TOP_KM`.

Resolution: `MULTIPLICATION_FACTOR` sets the lateral count (480·m control
points around the circle; air elements are 2πR/(240·m) wide). It is chosen
from `MIN_PERIOD_S` so air elements are ≤ 0.8·c_min·T (5 points per
wavelength, c_min ≈ 207 m/s). Interior radial counts are
`INTERIOR_RADIAL_FACTOR`·m (default 8, originally 56/24/40). This matters:
the layer shift squeezed 56·m radial elements into the 738 km upper
mantle, making elements 10× flatter than wide and setting the time step.
The air radial count is about twice the lateral spacing. The generator
repeats the mesher's own checks (even/multiple-of-4 counts) before
writing.

### Model table (`DATA/venus_model.txt`)

```text
# comments
N
region radius_m rho_kg_m3 vp_m_s vs_m_s qkappa qmu     (N rows)
```

Rows are grouped by region and sorted by radius within it. The Fortran
reader interpolates only between rows of an element's own region and
clamps at region ends, so it never blends across a discontinuity. vs = 0
marks fluid. Q = 9999 means no attenuation. Solid Q is forced to 9999
unless `SOLID_ATTENUATION = True`. Air rows are generated every 1 km:
density is interpolated in log space, sound speed linearly.

### Other choices

- **"No atmosphere" reference** = the same mesh with air density ×10⁻³.
  Verified: air pressure scales by 1.005×10⁻³. At a 60 s period, ground
  motion at 10° changes by only 0.5% between air and no air.
- **Time stepping**: Newmark (`time_stepping_scheme = 1`). `DT` = 0.8 ×
  the estimated maximum. The estimate is calibrated against the solver's
  "Maximum suggested time step" (`DT_CALIBRATION = 1.44`), and it matched
  the solver at m = 18 and m = 22 to within 1%.
- **Source**: a vertical force or an explosion on the axis (x = 0). These
  are the only valid axisymmetric sources; a horizontal force on the axis
  is invalid. Default depth is 10 km.
- **Stations**: ground stations are `G<deg>`, 20 m below the surface. Air
  stations are `A<deg>H<alt_km>` at fixed altitudes. Angles are measured
  from the source (+z axis).

## What was wrong with the original setup (for history)

The starting point was a modified copy of
`EXAMPLES/applications/axisymmetric_examples/global_Earth_ak135f`. Issues
found and fixed:

1. **Region flags.** They were renumbered in the mesher but not in the
   stock AK135F model routine, so the atmosphere became a dense solid.
   The student later wrote her own routine. It is now replaced by the
   table reader.
2. **Air sound speed** was a constant 1000 m/s. It is now a VIRA-like
   profile, 406 → 207 m/s.
3. **The top of the air was a rigid lid.** The free-surface list was
   empty, and an acoustic edge with no condition is rigid. PML was
   switched on with no CPML file, and SPECFEM2D's CPML is Cartesian, so it
   cannot follow a circular boundary. Stacey and rotated PML exit
   immediately with AXISYM. Together with the density drop, the air rang.
   It is now a Q sponge.
4. **The air was under-resolved**, with 80 km × 6.7 km elements. The
   lateral doubling sat directly beneath the solid surface.
5. **Mesher and table radii disagreed** by 1 m at the CMB; a "CMB FIX"
   clamp papered over it. Radii now come from the table.
6. **Solver settings**: LDDRK with DT = 0.01 and 10⁶ steps; receivers
   placed exactly on the fluid–solid interface; a horizontal force on the
   axis.

## SPECFEM2D bug: fluid attenuation in AXISYM (fixed by patch)

In `src/specfem2D/compute_forces_viscoacoustic.F90` (near lines 371–395,
the `ATTENUATION_VISCOACOUSTIC` AXISYM block), the memory-variable term is
added as `wzgll*wxgll*tempx3` with **no radius weight**. The stiffness and
mass terms both carry r (`coord(1,iglob)`, or `r_xiplus1(i,j)` on axis
elements), and so does `rmass_inverse_e1`. Fluid attenuation is therefore
about r ≈ 6×10⁶ times too weak, which is effectively off.

The patch multiplies the `tempx3` and `tempx3_e1` terms by
`coord(1,iglob)` off the axis and by `r_xiplus1(i,j)` on it. It applies
cleanly to SPECFEM2D `devel` at 47b2b79d. It has not been reported
upstream yet; the PI may want to open a SPECFEM2D issue.

**Verification** (45 s period, a column of stations 0.2° from the axis):
- Ratio of the sponge run to a no-sponge run with air extended to 200 km:
  - without the patch: 1.00 everywhere;
  - with the patch: 1.00 below 70 km, then 0.84 at 80 km, 0.55 at 90 km,
    and about 0.3 at 96–98 km.
- Echo off the top, measured at 1–65 km: −10 to −14 dB with the patch.

Caveats:
- The decay is about twice exp(−π f t*) and nearly independent of
  frequency. This is likely a limit of the three-SLS fit at Q ≈ 3.
- A residual echo of 10–20% of the direct wave remains at 30–65 km.
  Possible improvements: a thicker or gentler (e.g. quadratic) ramp.

Until a SPECFEM2D build has the patch, use `SPONGE = False` and analyze
each record only before the top echo returns. `run_summary.txt` gives the
one-way vertical air travel time, about 360 s, so the echo returns about
720 s after each arrival.

## Cost and performance

- Measured on an Apple-silicon laptop: about 1.0–1.3 µs per element per
  time step per core (Newmark, NGLL = 5), plus about 15% with fluid
  attenuation on.
- Cost scales as f_max³ (each halving of the period costs 8×).
- 20 s period, m = 50: 1.72 M elements, DT ≈ 0.012 s, 2000 s simulated →
  about 90–275 core-hours, roughly 0.5–1.6 h on 168 cores.
- The mesher at m = 50 took 17 s and 1.2 GB of RAM, and its element count
  matched the generator's prediction exactly.
- Laptop runs: use `MIN_PERIOD_S` of 60–100 s. Example: 100 s, m = 10,
  71k elements, took 82 s for 7,200 steps on 8 cores.

## Current state (20 Nov 2026)

- `venus_config.py` settings: `EXPERIMENT = "venus_T60"`,
  `MIN_PERIOD_S = 60`, `SPONGE = False`. The local SPECFEM bin path is a
  placeholder to fill in.
- The student is running long-period laptop tests while the cluster is
  busy.
- Verified end to end on a scratch build:
  - the generator;
  - both run scripts, including the 1-process local run;
  - the binary check;
  - the model reader (compiled with SPECFEM2D's strict flags:
    `-std=f2008 -pedantic`);
  - the attenuation patch.
- Not yet done:
  - production 20 s runs;
  - inner-core production runs (`venus_config_inner_core.py`);
  - reporting the AXISYM bug upstream.

## Inner-core layout (implemented 23 Nov 2026)

Implemented from the plan below, with these choices:

- `INNER_CORE_RADIUS_KM`: None = default layout (unchanged output,
  checked by diffing old and new generator runs); `"table"` = take the
  inner core from an interior table whose core has vs > 0 rows; a number
  = split an all-liquid core there and use constant
  `INNER_CORE_RHO/VP/VS`.
- Model: V5-Th from the PI (`data/models/V5-Th_ICB_2500.out`), CMB
  3183.3 km, ICB 2493.4 km, IC vp 11.0, vs 3.5 km/s. Its air rows are
  ignored (VIRA air kept: within ~4% of it at the ground).
- `CORE_CASES` × `CASES` → folders `<air>_<core>`. Liquid control:
  `INNER_CORE_LIQUID = "outer_core"` (OC value at the ICB, constant to the
  center) or `"same_bulk"` (same rho and K, mu = 0).
- `radius_cube` is patched in every mesher (650 km in the default layout,
  `650/1221 · R_ICB` in the inner-core layout).
- Radial counts in the inner-core layout come from
  `INTERIOR_ASPECT_RATIO` (radial/lateral element size, 1.5), not
  `INTERIOR_RADIAL_FACTOR`.
- Mesher constraint found: `x3vol` is dimensioned with
  `nspec_rad_ICB_to_doubling_OC` but also used for the zone above it, so
  `nspec_rad_doubling_OC_to_CMB` must be ≤ it (EQUIVALENCE hides the
  overrun). The generator enforces this.
- `estimate_dt` now takes, per zone, 0.7 × min(radial, lateral element
  size); the central cube term is unchanged, so default-layout DT is
  identical.
- `templates/Par_file` gained `save_ASCII_kernels`, required by the
  laptop SPECFEM2D (785f9d93); the student had added it by hand to runs.

## Original plan for the inner-core layout

No Fortran changes are needed. The mesher's original Earth topology
already has a spherical ICB with fluid outside and solid inside. Map the
zones differently:

| Mesher zone | Flag | New meaning |
|---|---|---|
| d670 → surface | 4 | Air (unchanged) |
| CMB → d670 | 1 | Whole mantle: surface → Venus CMB (3352.75 km) |
| ICB → doubling → CMB | 2 | Liquid outer core |
| Cube → ICB + central cube | 3 | Inner core: solid, or liquid as a control |

Implementation plan for `make_venus_run.py` and `venus_config.py`:

1. **New config settings.** Add to `DEFAULTS`:
   - `INNER_CORE_RADIUS_KM = None`: None keeps the current layout (whole
     core liquid, two mantle regions).
   - `INNER_CORE_SOLID = True`.
   - Inner-core properties, e.g. `INNER_CORE_RHO`, `INNER_CORE_VP`,
     `INNER_CORE_VS`, `INNER_CORE_QKAPPA`, `INNER_CORE_QMU`. These are
     scenario values, since Venus's inner core is unconstrained. Constant
     values or a small table both work.
   - Optionally a `CORE_CASES` dict, analogous to `CASES`, to make solid
     and liquid runs on the same mesh.
2. **Region assignment** when `INNER_CORE_RADIUS_KM` is set:
   - Region 1 = all mantle rows (merge current regions 1 and 2; the table
     is continuous at the 5313.568 km mid-mantle boundary).
   - Region 2 = core rows from the inner-core radius to the CMB (fluid).
   - Region 3 = inner-core rows (0 → inner-core radius).
   - Interpolate the core profile at the inner-core radius and duplicate
     that row on both sides.
   - For the liquid control, region 3 takes the outer-core properties
     with vs = 0.
3. **Mesher radii**: `R670` = surface, `RCMB` = CMB, `RICB` = inner-core
   radius. The outer-core doubling radius follows automatically
   (`RICB + 0.56·(RCMB − RICB)`). Re-check the radial counts: the mantle
   zone now spans about 2700 km, so `nspec_rad_CMB_670` may need a larger
   factor for well-shaped elements.
4. **`write_nummaterial`**: the fluid regions become `{2, 4}`, and region 3
   is fluid only for the liquid control. `FLUID_REGIONS` is currently
   hard-coded as `(REG_CORE, REG_AIR)`.
5. **`estimate_dt`**: update the zone-to-velocity mapping. The central
   cube and the shell around it are now inner core; a low inner-core vs
   does not affect dt, but the zone widths change.
6. **Constraint.** The central cube (`radius_cube = 650 km`, inflated
   outward) must fit inside the inner core. Require
   `INNER_CORE_RADIUS_KM ≥ ~1000`, or patch `radius_cube` in the mesher
   template through the generator.
7. **Test** at a low resolution (100–150 s) with a scratch SPECFEM2D:
   - the solver log shows the expected acoustic/elastic element counts;
   - the run is stable;
   - solid and liquid runs differ only at stations whose paths cross the
     inner core.

   Inner-core phases are long-period and deep, so run without air or with
   the no-air density scale.

## Pointers into SPECFEM2D (`devel`, 47b2b79d)

- `src/specfem2D/read_external_model.f90`: calls `define_external_model`,
  then checks that vs = 0 agrees with the Nummaterial domain. That check
  produces the error "acoustic material … redefined as non acoustic".
- `src/specfem2D/axisymmetric_routines.f90:103-106`: rejects rotated PML
  and Stacey with AXISYM.
- `src/specfem2D/pml_init.F90:758`: CPML damping from a Cartesian bounding
  box.
- `src/specfem2D/prepare_attenuation.f90`: the absorption band is f0/10 to
  10·f0, because `USE_FIXED_ATTENUATION_ABSORPTION_BAND = .true.` in
  `setup/constants.h`.
- `src/shared/read_parameter_file.F90`: missing Par_file keys stop the
  run, so the Par_file template must match the SPECFEM2D version.
- The stock example
  `EXAMPLES/applications/axisymmetric_examples/axisymmetric_case_AXISYM_option`
  (internal mesher, AXISYM, a fluid layer, PML on three flat sides) is the
  template for a flat axisymmetric "chunk" alternative. That option was
  considered and set aside to keep the global geometry.
