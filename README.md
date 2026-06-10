# Seismic velocity conversion and tomographic filtering

This repository post-processes mantle convection simulation output into
synthetic seismic observables that can be compared directly — and fairly —
with global tomographic models. Given the PVTU output of a Firedrake/G-ADOPT
simulation, it:

1. converts temperature to absolute seismic velocities (Vs, Vp) using a
   thermodynamic look-up table, with gradient regularisation and an anelastic
   correction;
2. filters the Vs field through the resolution operators of **S12RTS**,
   **S20RTS**, and **S40RTS** (Ritsema et al. 1999, 2004, 2011);
3. filters the Vs and Vp fields through the resolution matrix of
   **LLNL-G3D-JPS** (Simmons et al. 2015, 2019);
4. resamples everything onto a regular lon/lat/depth grid as NetCDF;
5. renames and reshapes the variables into an analysis-ready convention.

On this branch (`mass-convert`) the five steps are orchestrated by a single
driver script, `pipeline.sh`, designed to batch-process many simulation
snapshots on a local machine. Each numbered step can also be run on its own.

```
Cratons_XMa.tar.gz  (one per time slice)
        |
        v  pipeline.sh ── extracts each tarball, then per snapshot:
XMa/output/output_0.pvtu
        |
        v  01_convert.sh        (convert_to_v.py, uses gdrift)
XMa/converted.vtu                Vs, Vp + dlnVs, dlnVp at every mesh node
        |
        +─ 02_srts_filter.sh    (srts_filter.py, uses srts)
        |  XMa/converted_srts_filtered.vtu      S12/S20/S40RTS-filtered Vs
        |
        +─ 03_llnl_filter.sh    (llnl_filter.py, uses llnltofi)
        |  XMa/converted_llnl_filtered.vtu      LLNL-G3D-JPS-filtered Vs, Vp
        |
        v  04_interpolate.sh    (ginterp)
XMa/converted{,_srts_filtered,_llnl_filtered}.nc   regular 360x181x129 grid
        |
        v  05_rename.sh         (rename.py)
<output_dir>/Cratons_XMa.nc
<output_dir>/Cratons_XMa_S40RTS_ToFi.nc
<output_dir>/Cratons_XMa_LLNL_ToFi.nc
```

The science behind each step is explained below, after the practical
instructions.

---

## Installation

The three core packages are listed in `requirements.txt`:

```bash
pip install -r requirements.txt
```

This installs:

- **srts** — <https://github.com/g-adopt/srts> — S12/S20/S40RTS tomographic
  filtering (from PyPI);
- **gdrift** — <https://github.com/g-adopt/g-drift> — thermodynamic
  conversion tables, gradient regularisation, and anelastic corrections
  (from GitHub);
- **llnltofi** — <https://github.com/g-adopt/llnltofi> — LLNL-G3D-JPS
  resolution-matrix filtering (from GitHub);
- **ginterp** — <https://github.com/g-adopt/g-interp> — mesh-to-grid
  interpolation for step 4 (from GitHub):

  ```bash
  pip install git+https://github.com/g-adopt/g-interp.git
  ```

The scripts additionally use `pyvista` (VTK input/output), `xarray` (NetCDF
handling), and `numpy`/`scipy`; make sure these are available in your
environment. `pipeline.sh`
calls GNU tar as `gtar` — on macOS install it with `brew install gnu-tar`;
on Linux, alias or symlink `gtar` to `tar`.

The model data files required by `srts`, `gdrift`, and `llnltofi` (look-up
tables, resolution operators, grid geometries — several GB in total) are
downloaded automatically and cached on first use, so the first run needs
internet access.

---

## Running the pipeline

### Batch mode (the usual way)

Point `pipeline.sh` at a directory of `Cratons_XMa.tar.gz` tarballs (one per
reconstruction age `X`, in Ma) and at a destination for the final NetCDF
files:

```bash
./pipeline.sh <base_dir> <output_dir> [scripts_dir] [--min-age N] [--max-age N]
```

- `base_dir` — directory containing the `Cratons_XMa.tar.gz` files; an
  `XMa/` working subdirectory is created next to each tarball;
- `output_dir` — destination for the final renamed `.nc` files;
- `scripts_dir` — location of the numbered scripts and Python helpers
  (defaults to the directory containing `pipeline.sh`);
- `--min-age` / `--max-age` — optionally restrict processing to tarballs
  within an age range (inclusive, in Ma).

Example:

```bash
./pipeline.sh /Volumes/Grey/firedrake_simulations/Cratons \
              ~/OneDrive/phd/firedrake-models/Cratons \
              --min-age 10 --max-age 50
```

For each tarball the driver extracts the archive, runs steps 1–5, and
reports a processed/skipped tally at the end. All intermediate artefacts
(`.vtu` and unrenamed `.nc`) are kept in the per-age working directory; only
the renamed NetCDF files land in `output_dir`.

### Running a single snapshot

Each step is a standalone script taking the working directory as its first
argument, so a single snapshot can be (re)processed step by step:

```bash
bash 01_convert.sh      <work_dir>
bash 02_srts_filter.sh  <work_dir>
bash 03_llnl_filter.sh  <work_dir>
bash 04_interpolate.sh  <work_dir>
bash 05_rename.sh       <work_dir> <output_dir>
```

`<work_dir>` must contain the simulation output at
`output/output_0.pvtu`. Steps 2 and 3 are independent of each other and can
run in either order (or concurrently); both require step 1, and step 4
requires all three VTU files.

---

## The physics, step by step

### 1. From temperature to seismic velocity (`convert_to_v.py`)

The simulation evolves a non-dimensional temperature field. To compare it
with seismology we must predict what seismic waves would measure in that
mantle: Vs and Vp as a function of temperature and depth. The script
dimensionalises coordinates and temperature (surface temperature 300 K,
temperature drop 3700 K, mantle depth 2891 km, non-dimensional outer radius
2.208 — constants at the top of the file), reading the
`FullTemperature_CG` and `Temperature_Deviation_CG` fields from the PVTU.

**Thermodynamic model.** Conversion uses the SLB_24 dataset (Stixrude &
Lithgow-Bertelloni 2024) with a pyrolite CFMASNaCr
(CaO–FeO–MgO–Al₂O₃–SiO₂–Na₂O–Cr₂O₃) bulk composition, accessed through
`gdrift.ThermodynamicModel`. This is a pre-computed look-up table: for each
(temperature, depth) pair it gives Vs, Vp, and density for the stable
mineral assemblage at those conditions, covering the full mantle from 0 to
2891 km.

**Gradient regularisation.** Phase transitions (olivine → wadsleyite →
ringwoodite → post-spinel, etc.) produce sharp velocity jumps in the raw
table as a function of temperature at fixed depth, so a naive conversion
produces artefacts wherever the temperature field crosses a phase boundary.
`gdrift.regularise_thermodynamic_table` anchors the conversion to the
spherically averaged temperature profile T̄(z) of the simulation itself —
computed here by averaging the mesh temperatures in 200 evenly spaced depth
bins — and maps temperature anomalies to velocity anomalies linearly about
the reference velocity V(T̄(z), z), with the velocity–temperature gradient
clipped to a physically plausible range. The result is smooth and free of
phase-transition artefacts.

**Anelastic correction.** The SLB_24 table gives *elastic* (infinite-
frequency) velocities, but real seismic waves travel at finite frequency,
and at high mantle temperatures anelastic attenuation causes dispersion:
the measured velocity is lower than the elastic one. The correction follows
Cammarano et al. (2003), implemented in `gdrift.CammaranoAnelasticityModel`
with the **Q6** quality-factor profile, and is largest (several per cent) in
the hot deep mantle where Q is low.

**Outputs.** The script writes `converted.vtu` containing three velocity
variants at every mesh node, so that the effect of each modelling choice can
be isolated downstream:

| Field(s)                   | Meaning                                            |
|----------------------------|----------------------------------------------------|
| `T`, `dT`                  | temperature (K) and temperature deviation (K)      |
| `Vs_unreg`, `Vp_unreg`     | raw table conversion (no regularisation)           |
| `Vs_lin`, `Vp_lin`         | gradient-regularised (linearised), elastic         |
| `Vs`, `Vp`                 | gradient-regularised **and** anelastically corrected |
| `dlnVs_*`, `dlnVp_*`       | relative anomalies (%) of each variant             |

The relative anomalies are computed by `_layer_mean.py` as the deviation
from the unweighted nodal mean within each depth layer of the extruded mesh.
Because the horizontal mesh is quasi-uniform (near-equal area per node), the
unweighted layer mean is a faithful approximation to the true spherical
mean, so no cos(latitude) weighting is needed at this stage.

### 2. S-RTS tomographic filtering (`srts_filter.py`)

A direct comparison between the synthetic Vs field and a tomographic model
is unfair: the tomographic model sees a blurred and damped version of the
Earth, determined by the source–receiver geometry and the inversion
regularisation. Tomographic filtering imposes the same blurring on the
synthetic field so that like is compared with like.

The S-RTS family parameterises Vs anomalies as spherical harmonics
horizontally and 21 vertical splines. Using the `srts` package, the script:

1. interpolates the mesh onto a regular 181 × 360 lat/lon grid at each
   depth layer by inverse-distance weighting (k = 4 neighbours, power 2);
2. expands each layer in spherical harmonics to degree 40
   (`SphericalHarmonicExpansion`);
3. projects the coefficient profiles onto the 21-knot spline depth basis
   (`DepthParameterization`) — a change of vertical basis, not yet a
   filtering step;
4. applies the resolution operator of each of S40RTS, S20RTS, and S12RTS,
   which truncates horizontal structure beyond each model's resolution
   (ℓmax = 40, 20, 12) and damps vertical structure the inversion cannot
   constrain;
5. synthesises the filtered coefficients back onto the regular grid and
   interpolates back to the mesh nodes layer by layer, restoring the layer
   mean so absolute velocities are preserved.

Both IDW weight matrices are precomputed once from the shared horizontal
structure of the extruded mesh, and the spherical harmonic operators are
reused across layers, so the whole filter runs in minutes.

**Outputs** (added to `converted_srts_filtered.vtu`): `Vs_reparam` — Vs
after the round trip through the spectral–spline basis but *without*
filtering (the correct unfiltered reference for filtered-vs-unfiltered
comparisons); `Vs_tofi` — S40RTS-filtered; `Vs_tofi_S20RTS` and
`Vs_tofi_S12RTS`; and the corresponding `dlnVs_*` percentage anomalies.

### 3. LLNL-G3D-JPS tomographic filtering (`llnl_filter.py`)

LLNL-G3D-JPS (Simmons et al. 2015, 2019) is a joint P- and S-wave model on
an irregular grid whose point density follows ray-path coverage. Its
resolution matrix **R** is the explicit least-squares operator

    R = (GᵀC_d⁻¹G + C_m⁻¹)⁻¹ GᵀC_d⁻¹G

where G is the sensitivity matrix, C_d the data covariance, and C_m the
model covariance (regularisation). Applying **R** to a synthetic slowness
anomaly gives what the LLNL inversion would have recovered had the Earth
looked like the simulation. Using the `llnltofi` package, the script:

1. projects Vs and Vp from the mesh onto the ~10⁶-point LLNL grid with a
   layer-aware IDW (`project_onto_grid`);
2. converts to slowness and forms the anomaly δs = s − s₁D relative to the
   layer-mean 1-D reference;
3. applies the resolution matrix as a single sparse matrix–vector product;
4. converts back to velocity and back-projects to the mesh nodes
   (`project_from_grid`), mirroring the forward interpolation on the same
   layer geometry.

No separate amplitude treatment is needed for Vs versus Vp: **R** acts on
slowness anomalies irrespective of wave type.

**Outputs** (added to `converted_llnl_filtered.vtu`): `Vs_reparam`,
`Vp_reparam` — the unfiltered round trip through the LLNL grid (the correct
unfiltered reference for this filter); `Vs_tofi`, `Vp_tofi` — filtered; and
the corresponding `dlnVs_*`/`dlnVp_*` percentage anomalies.

### 4. Interpolation to a regular grid (`04_interpolate.sh`)

The VTU files live on the unstructured finite-element mesh, which is
convenient for computation but awkward for analysis. Step 4 uses
[`ginterp`](https://github.com/g-adopt/g-interp) to resample all three onto
a regular 360 × 181 × 129
(longitude × latitude × depth) spherical grid spanning non-dimensional
radii 1.208–2.208 (CMB to surface), writing one NetCDF file per VTU.

### 5. Renaming for analysis (`05_rename.sh`, `rename.py`)

The final step subsets and renames the NetCDF variables into the convention
used by the downstream plotting and analysis workflows, and tidies the
coordinates: dimensions are reordered to (r, lat, lon), longitudes are
wrapped to [−180°, 180°) and sorted, the radial coordinate is
dimensionalised to metres, and a `depth` coordinate (km) is added. The
variable mapping is:

| Input variable  | Output variable          | Meaning                                  |
|-----------------|--------------------------|------------------------------------------|
| `T`, `dT`       | `T`, `dT`                | temperature, temperature deviation       |
| `dlnVs_lin`     | `dlnVs_lin_percent`      | regularised, elastic                     |
| `dlnVs`         | `dlnVs_linan_percent`    | regularised + anelastic                  |
| `dlnVs_reparam` | `dlnVs_reparam_percent`  | unfiltered reference (reparameterised)   |
| `dlnVs_tofi`    | `dlnVs_tofi_percent`     | tomographically filtered                 |

(and likewise for `dlnVp_*` where present). Three files are written per
snapshot, named after the working directory: `<TAG>.nc` (conversion only),
`<TAG>_S40RTS_ToFi.nc`, and `<TAG>_LLNL_ToFi.nc`. These are the primary
data products for making maps, radial profiles, and power spectra.

---

## Legacy HPC scripts

Two scripts from the Gadi (NCI) workflow that preceded this branch are kept
for staging and archiving data on the HPC system; they are not called by
`pipeline.sh`:

- `00_stage.sh` — PBS `copyq` job that extracts a simulation tarball into
  the scratch working directory and sanity-checks the expected PVTU;
- `archive.sh` — PBS `copyq` job that tars a run's six pipeline artefacts
  into the project archive on `/g/data`.

---

## Key references

- Stixrude & Lithgow-Bertelloni (2024) — SLB_24 thermodynamic dataset
- Cammarano et al. (2003) — anelastic velocity corrections and Q profiles
- Ritsema et al. (1999, 2004, 2011) — S12RTS, S20RTS, S40RTS
- Simmons et al. (2015, 2019) — LLNL-G3D-JPS and its resolution matrix
- Ritsema et al. (2007) — tomographic filtering of geodynamic models