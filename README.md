# Seismic velocity conversion and tomographic filtering

This repository contains the scripts needed to post-process mantle convection
simulation output into synthetic seismic observables that can be directly
compared with global tomographic models. The pipeline has five steps---four
science steps plus one step to harmonise variable names---which are performed
automatically on tarballs with `pipeline.sh`.

```
simulation output tarball (.tar.gz of *.pvtu pieces)
        |
        v  00_stage.sh                       (one-off untar onto /scratch)
${NAME}_output/output/output_0.pvtu
        |
        v  01_convert.sh
${NAME}_converted.vtu                        Vs, Vp, dlnVs, dlnVp at every mesh node
        |
        +--> 02_srts_filter.sh
        |    ${NAME}_converted_srts_filtered.vtu     S40/S20/S12 RTS filtered Vs + dlnVs
        |
        +--> 03_tofi_filter.sh
        |    ${NAME}_converted_tofi_filtered.vtu     LLNL-G3D-JPS filtered Vs, Vp + dlnVs, dlnVp
        |
        +--> 04_interpolate.sh
             ${NAME}_converted*.nc                   all three on a regular lon/lat/depth grid
```

Every script is generic on a single run identifier `NAME` passed via
`qsub -v NAME=...`, so the same pipeline runs against any number of
simulations side by side without editing source.  Active run names and the
tarballs they came from are listed in `RUNS.md`.

The science behind each step is explained below, followed by practical
instructions for running the jobs.

---

## Installation

The three Python packages this pipeline depends on are listed in
`requirements.txt`.  Install them into your environment with:

```bash
pip install -r requirements.txt
```

This pulls `srts` from PyPI and installs `gdrift` and `llnltofi` directly from
their GitHub repositories:

- **srts** - PyPI package for S-RTS tomographic filtering
- **gdrift** - <https://github.com/g-adopt/g-drift> - thermodynamic conversion tables and anelastic corrections
- **llnltofi** - <https://github.com/g-adopt/llnltofi> - LLNL-G3D-JPS resolution matrix filtering

On Gadi the packages are already installed under the `xd2` project's shared
Python path, so you do not need to re-install them - just load the module and
run the PBS scripts as described in the usage section below.

---

## The physics, step by step

### 1. From temperature to seismic velocity

The simulation evolves a non-dimensional temperature field through the mantle.
To compare it with seismology we need to predict what seismic waves would
measure in that mantle, that is, we need Vs (shear-wave velocity) and Vp
(compressional-wave velocity) as a function of temperature and pressure.

**Thermodynamic model.** We use the SLB_24 dataset (Stixrude &
Lithgow-Bertelloni 2024) with a pyrolite CFMASNaCr (CaO-FeO-MgO-Al2O3-SiO2-
Na2O-Cr2O3) bulk composition.  This is a pre-computed thermodynamic look-up table:
for each (temperature, depth) pair it gives Vs, Vp, and density, derived from
mineral-physics equations of state for the stable phase assemblage at those
conditions.  The table covers the full mantle from 0 to 2891 km. These are
coefficients coming from the above study and then thrown into a numerical code
that optimises different parameters for mineral assemblage to find what the most
stable mineralogy is. Then using theoretical expectations we compute all the
thermodynamic parameters.

**Why regularisation?**  Phase transitions (olivine -> wadsleyite -> ringwoodite
-> post-spinel, etc.) produce sharp velocity jumps in the raw table as a
function of temperature at a fixed depth.  In a convecting mantle the average
temperature at any depth is not zero; it follows a geotherm, so a naive
conversion produces artefacts wherever the laterally-averaged temperature
crosses a phase boundary.  `regularise_thermodynamic_table` anchors the
conversion to the actual spherically-averaged temperature profile extracted from
the simulation mesh.  Concretely, it computes $\bar{T}(z)$ from the simulation,
evaluates the reference velocity $V_s^{\rm ref}(\bar{T}(z), z)$, and then maps
temperature anomalies $\delta T$ to velocity anomalies $\delta V_s$ linearly
around that reference.  The result is smooth and free of phase-transition
artefacts.

**Anelastic correction.** The SLB_24 table gives *elastic* velocities, i.e.
what you would measure at infinite frequency.  Real seismic waves travel at
roughly 1 Hz, and at the high temperatures in the deep mantle this matters:
anelastic attenuation causes velocity dispersion, meaning the actual seismic
velocity is measurably lower than the elastic one.  The correction follows
Cammarano et al. (2003) and uses the Q3 quality-factor profile, which is
calibrated against observed 1-D seismological reference models.  The
corrected velocity is

$$
V_s^{\rm anel}(T, z) \approx V_s^{\rm el}(T, z)
\left[1 - \frac{\cot(\alpha\pi/2)}{2 Q(z)}\right]
$$

where $\alpha \approx 0.26$ is the frequency exponent and $Q(z)$ is the depth-
dependent quality factor.  The correction is largest (several percent) in the
deep lower mantle where temperatures are high and Q is low.

**Input / output.**  The script reads one PVTU file (the parallel VTK format
written by Firedrake/G-ADOPT), dimensionalises coordinates and temperature using
the constants at the top of `convert_to_vs.py`, and writes a single VTU with
fields `Temperature_K`, `Vs`, `Vp`, `dlnVs`, and `dlnVp` at every mesh node.

The `dlnVs` and `dlnVp` fields are the linearised seismological perturbations
relative to the depth-mean velocity,

$$
\delta\ln V \equiv \frac{V - \langle V\rangle(z)}{\langle V\rangle(z)} \times 100\%,
$$

where $\langle V\rangle(z)$ is the spherical-shell average at the depth of each
node.  Mesh nodes lie on discrete radial layers, and on the quasi-uniform
icosahedral horizontal mesh each node represents nearly equal area, so the
unweighted nodal mean within a layer is already a faithful area-weighted
spherical mean — no $\cos\phi$ weighting is needed at this stage.  The fields
are carried through steps 2–4 unchanged, so the regularly gridded NetCDF
output of step 4 inherits the same $\delta\ln V$ values without any additional
averaging logic.

---

### 2. S-RTS tomographic filtering (srts_filter.py)

A direct comparison between the raw synthetic Vs field and a global tomographic
model is "unfair" because the tomographic model does not see the true Earth. It
sees a blurred and damped version of it, determined by the geographic
distribution of seismic sources and receivers and by the choice of inversion
regularisation.  Tomographic filtering replicates that effect on the synthetic
field so that like is compared with like.

The SRTS family - S12RTS (Ritsema et al. 1999), S20RTS (Ritsema et al. 2004),
S40RTS (Ritsema et al. 2011) - parameterises Vs anomalies as a sum of spherical
harmonics horizontally and 21 splines vertically.  The filtering procedure is:

1. **Mesh -> regular grid (IDW).** The unstructured simulation mesh is
   interpolated onto a 181 x 360 (latitude x longitude) regular grid at each
   depth layer using inverse-distance weighting.

2. **Spherical harmonic expansion.** The velocity anomaly on the regular grid
   is expanded in spherical harmonics up to degree $\ell_{\max}=40$.  This
   gives a set of spectral coefficients $c_{\ell m}(z)$ at each depth.

3. **Depth reparameterisation.** The continuous depth profile of coefficients
   is projected onto the 21 B-spline basis functions used by each S-RTS model.
   This is a change of vertical basis, not a filtering step.

4. **Resolution filter.** Each S-RTS model provides its own resolution operator
   in spectral-spline space.  It encodes both the maximum resolved wavelength
   (given by $\ell_{\max}$: 12, 20, or 40) and the vertical smoothing imposed
   by the inversion.  Applying the filter truncates horizontal structure beyond
   the model's resolution and damps vertical structure that the inversion could
   not constrain.

5. **Synthesis + grid -> mesh (layered).** The filtered coefficients are
   synthesised back to the regular grid and then interpolated back to the mesh
   nodes layer-by-layer.  The depth-layer mean is restored before the final
   output so that absolute velocities (not just anomalies) are preserved.

The output VTU adds three filtered fields alongside the original Vs:
`Vs_S40RTS`, `Vs_S20RTS`, and `Vs_S12RTS`, plus their linearised
perturbations `dlnVs_S40RTS`, `dlnVs_S20RTS`, and `dlnVs_S12RTS`.  Each
`dlnVs_*` is referenced to the depth-mean of its own filtered field,
following the standard tomographic-anomaly convention.

---

### 3. LLNL-G3D-JPS tomographic filtering (tofi_filter.py)

LLNL-G3D-JPS (Simmons et al. 2012, 2019) is a joint P- and S-wave tomographic
model parameterised on an irregular grid with a geographic point distribution
that follows ray-path density.  Its resolution matrix $\mathbf{R}$ is the
explicit least-squares solution to

$$
\mathbf{R} = (\mathbf{G}^T \mathbf{C}_d^{-1} \mathbf{G}
+ \mathbf{C}_m^{-1})^{-1} \mathbf{G}^T \mathbf{C}_d^{-1} \mathbf{G}
$$

where $\mathbf{G}$ is the ray-path sensitivity matrix, $\mathbf{C}_d$ the data
covariance, and $\mathbf{C}_m$ the model covariance (regularisation).  Applying
$\mathbf{R}$ to a synthetic slowness vector gives what the LLNL inversion would
have recovered had the Earth looked like the simulation.

The filtering steps are:

1. **Mesh -> LLNL grid (IDW).** A layer-by-layer IDW maps the simulation Vs and
   Vp onto the geographic grid points of the LLNL model.  The LLNL grid
   distinguishes upper-mantle/transition-zone layers from lower-mantle layers,
   each with different horizontal resolutions.

2. **Convert to slowness anomaly.** The depth-dependent 1D reference is
   estimated from the layer-mean slowness on the LLNL grid, and the anomaly
   $\delta s = s - s_{1\rm D}$ is formed.

3. **Apply resolution matrix.** $\delta s_{\rm filtered} = \mathbf{R}\,\delta s$.
   No separate amplitude scaling is applied for Vs vs Vp because $\mathbf{R}$
   acts on slowness anomalies irrespective of wave type.

4. **Recover velocity + back-projection.** The filtered slowness anomaly is
   added back to the 1D reference and converted to velocity.  A layered
   back-projection (`llnltofi.interpolation.project_from_grid`) maps the
   result back to the simulation mesh nodes layer-by-layer, mirroring the
   forward IDW step on the same LLNL layer geometry.

The output VTU adds `Vs_filtered` and `Vp_filtered`, plus the linearised
perturbations `dlnVs_filtered` and `dlnVp_filtered`, each referenced to the
depth-mean of the filtered field itself.

---

### 4. Interpolation to a regular grid (04_interpolate.sh)

The filtered VTU files live on the unstructured finite-element mesh, which is
convenient for computation but awkward for analysis and plotting.  Step 4 uses
`ginterp` to resample all three outputs onto a 360 x 181 x 129
(longitude x latitude x depth) regular grid and writes them as NetCDF files.
These are what you load for making maps, radial profiles, and power spectra.

---

## Running the pipeline on Gadi

### Prerequisites

The required Python packages are available through the Firedrake module on
Gadi.  No additional installation is needed beyond what is already set up in
the `xd2` project — the PBS scripts `module load firedrake/...` and prepend
the project-shared `gdrift`, `srts`, `llnltofi`, and `ginterp` paths to
`PYTHONPATH` for you.

### Setup

Clone this repository somewhere under your scratch space:

```bash
cd /scratch/xd2/USERNAME
git clone <repo-url> kat-conversion
cd kat-conversion
```

Replace `USERNAME` with your Gadi username throughout.

### How the scripts are parametrised

Every step is generic on a run identifier `NAME` (and, for staging, an
`INPUT_TAR` path), passed via `qsub -v`.  Files for a given run all share
the prefix `${NAME}_`, so multiple simulations coexist in the same working
directory without collision.  The naming convention is:

| Step              | Reads                                                    | Writes                                      |
|-------------------|----------------------------------------------------------|---------------------------------------------|
| `00_stage.sh`     | `INPUT_TAR` (tarball of pvtu pieces)                     | `${NAME}_output/output/output_0.pvtu`       |
| `01_convert.sh`   | `${NAME}_output/output/output_0.pvtu`                    | `${NAME}_converted.vtu`                     |
| `02_srts_filter.sh` | `${NAME}_converted.vtu`                                | `${NAME}_converted_srts_filtered.vtu`       |
| `03_tofi_filter.sh` | `${NAME}_converted.vtu`                                | `${NAME}_converted_tofi_filtered.vtu`       |
| `04_interpolate.sh` | the three VTUs above                                   | `${NAME}_converted{,_srts,_tofi}_filtered.nc` |

`01_convert.sh` also accepts an `INPUT_PVTU=...` override if you want to
point at a pvtu that wasn't placed by `00_stage.sh`.

### End-to-end submission

Pick a `NAME` and the source tarball, then submit the five jobs with PBS
dependencies so they queue immediately and start in order as soon as each
predecessor finishes.  `02_srts_filter.sh` and `03_tofi_filter.sh` both
depend on `01_convert.sh` and run in parallel.

```bash
NAME=C39_3e22_MuT
TAR=/g/data/xd2/sg8812/kat-conversion-archive/0Ma_C39_3e22_MuT_Output.tar.gz

STAGE=$(qsub  -v NAME=$NAME,INPUT_TAR=$TAR                                 00_stage.sh)
CONV=$(qsub   -v NAME=$NAME -W depend=afterok:$STAGE                       01_convert.sh)
SRTS=$(qsub   -v NAME=$NAME -W depend=afterok:$CONV                        02_srts_filter.sh)
TOFI=$(qsub   -v NAME=$NAME -W depend=afterok:$CONV                        03_tofi_filter.sh)
INTERP=$(qsub -v NAME=$NAME -W depend=afterok:$SRTS:$TOFI                  04_interpolate.sh)
```

If the input has already been staged (i.e. `${NAME}_output/` exists), skip
the `00_stage.sh` line and drop `-W depend=afterok:$STAGE` from the
`01_convert.sh` line.

### Resource budgets

| Step | Queue   | ncpus | mem    | walltime |
|------|---------|------:|-------:|---------:|
| 00   | copyq   | 1     | 8 GB   | 4 h      |
| 01   | normal  | 1     | 128 GB | 4 h      |
| 02   | normal  | 1     | 64 GB  | 2 h      |
| 03   | normal  | 1     | 64 GB  | 2 h      |
| 04   | normal  | 1     | 128 GB | 4 h      |

### Checking job status

```bash
qstat -u USERNAME          # all your jobs
qcat -o <jobid>            # stdout of a completed job
qstat -fx <jobid>          # full record incl. exit status
```

### Tracking and archiving runs

`RUNS.md` records each `NAME` together with its source tarball and any
relevant notes (viscosity, rheology, model series).  Add a row whenever you
launch a new run.  Once `04_interpolate.sh` for a given `NAME` finishes,
copy the six `${NAME}_converted*` files (three `.vtu` plus three `.nc`)
into `/g/data/xd2/sg8812/kat-conversion-archive/` as
`kat_${NAME}_artefacts.tar.gz` from a `copyq` job so `/scratch` can be
reclaimed.

---

## Key references

- Stixrude & Lithgow-Bertelloni (2005, 2024) - SLB thermodynamic framework
- Cammarano et al. (2003) - Anelastic velocity corrections, Q3 profile
- Ritsema et al. (1999, 2004, 2011) - S12RTS, S20RTS, S40RTS
- Simmons et al. (2012, 2019) - LLNL-G3D-JPS and resolution matrix
