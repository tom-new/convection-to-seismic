"""
tofi_filter.py — Tomographic filtering of Vs and Vp using the
LLNL-G3D-JPS resolution matrix with layered interpolation.

Usage:
    python tofi_filter.py <input.vtu> <output.vtu>

The input is the output of convert_to_vs.py: a VTU (or PVTU) with Vs and Vp
as point-data fields (m/s) and non-dimensional Cartesian coordinates
(RMAX=2.208, D_KM=2891).  The script adds Vs_filtered and Vp_filtered
fields and writes a new VTU.

Forward projection (mesh -> LLNL grid) uses the layer-by-layer IDW in
llnltofi.interpolation.project_onto_grid.  Back-projection (LLNL -> mesh)
uses the layered interpolation in llnltofi.interpolation.project_from_grid.
"""

import sys
from pathlib import Path

import numpy as np
import pyvista as pv

import llnltofi
from llnltofi.interpolation import project_onto_grid, project_from_grid
from llnltofi._constants import (
    N_LAYERS_UM_TZ,
    N_LAYERS_LM,
    N_POINTS_UM_TZ,
    N_POINTS_LM,
    R_EARTH_KM,
)

from _layer_mean import dln_percent_by_layer

RMAX = 2.208
D_KM = 2891.0


# -- coordinate helpers --------------------------------------------------------

def nondim_to_spherical(coords):
    """Non-dimensional Cartesian -> (gc_lat_deg, lon_deg, radius_km).

    This is the format expected by project_onto_grid and project_from_grid.
    """
    r = np.linalg.norm(coords, axis=1)
    depth_km = (RMAX - r) * D_KM
    radius_km = R_EARTH_KM - depth_km
    unit = coords / r[:, np.newaxis]
    lat_deg = np.degrees(np.arcsin(np.clip(unit[:, 2], -1.0, 1.0)))
    lon_deg = np.degrees(np.arctan2(unit[:, 1], unit[:, 0]))
    return np.column_stack([lat_deg, lon_deg, radius_km])


# -- filtering helpers ---------------------------------------------------------

def layer_mean_1d(slowness):
    """Depth-dependent 1D reference from layer-wise mean slowness."""
    n_um_tz = N_LAYERS_UM_TZ * N_POINTS_UM_TZ
    s_um = slowness[:n_um_tz].reshape(N_LAYERS_UM_TZ, N_POINTS_UM_TZ)
    s_lm = slowness[n_um_tz:].reshape(N_LAYERS_LM, N_POINTS_LM)
    return np.concatenate([
        np.repeat(s_um.mean(axis=1), N_POINTS_UM_TZ),
        np.repeat(s_lm.mean(axis=1), N_POINTS_LM),
    ])


def apply_filter(velocity_on_llnl, model):
    """Filter one velocity field through the resolution matrix.

    Converts to slowness, removes the layer-mean 1D reference, applies R,
    then recovers velocity.
    """
    s = 1.0 / velocity_on_llnl
    s_1d = layer_mean_1d(s)
    model.values = s - s_1d
    ds_filtered = model.apply()
    return 1.0 / (s_1d + ds_filtered)


# -- main ----------------------------------------------------------------------

def main():
    if len(sys.argv) < 3:
        print(f"Usage: {sys.argv[0]} <input.vtu> <output.vtu>")
        sys.exit(1)

    input_path = sys.argv[1]
    output_path = sys.argv[2]

    # -- 1. Load LLNL model ---------------------------------------------------
    print("Loading LLNL resolution model...")
    model = llnltofi.ResolutionModel()
    _ = model.R                              # ensure R.npz is present
    print(f"  {model.n_model:,} grid points, {model.n_layers} layers")

    # -- 2. Load simulation mesh -----------------------------------------------
    print(f"\nReading {input_path} ...")
    mesh = pv.read(input_path)
    coords = np.asarray(mesh.points)
    vs = np.asarray(mesh.point_data["Vs"])
    vp = np.asarray(mesh.point_data["Vp"])
    print(f"  {len(vs):,} mesh points")

    # -- 3. Convert mesh coordinates -------------------------------------------
    print("Converting coordinates...")
    sph_coords = nondim_to_spherical(coords)

    # -- 4. Forward IDW: mesh -> LLNL grid -------------------------------------
    print("\nForward IDW: Vs mesh -> LLNL grid...")
    vs_on_llnl = project_onto_grid(sph_coords, vs, model)
    print("Forward IDW: Vp mesh -> LLNL grid...")
    vp_on_llnl = project_onto_grid(sph_coords, vp, model)

    # -- 5. Apply resolution filter --------------------------------------------
    print("\nFiltering Vs through resolution matrix...")
    vs_filtered_llnl = apply_filter(vs_on_llnl, model)
    print("Filtering Vp through resolution matrix...")
    vp_filtered_llnl = apply_filter(vp_on_llnl, model)

    # -- 6. Layered back-projection: LLNL grid -> mesh -------------------------
    print("\nBack-projection (layered): Vs LLNL -> mesh...")
    vs_filtered_mesh = project_from_grid(vs_filtered_llnl, sph_coords, model)
    print("Back-projection (layered): Vp LLNL -> mesh...")
    vp_filtered_mesh = project_from_grid(vp_filtered_llnl, sph_coords, model)

    # -- 7. Linearised seismological perturbation -----------------------------
    print("\nComputing dlnVs_filtered and dlnVp_filtered ...")
    depth_km = (RMAX - np.linalg.norm(coords, axis=1)) * D_KM
    dlnvs_filtered = dln_percent_by_layer(vs_filtered_mesh, depth_km)
    dlnvp_filtered = dln_percent_by_layer(vp_filtered_mesh, depth_km)

    # -- 8. Write output -------------------------------------------------------
    mesh.point_data["Vs_filtered"] = vs_filtered_mesh
    mesh.point_data["Vp_filtered"] = vp_filtered_mesh
    mesh.point_data["dlnVs_filtered"] = dlnvs_filtered
    mesh.point_data["dlnVp_filtered"] = dlnvp_filtered

    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    print(f"\nWriting {output_path} ...")
    mesh.save(output_path)

    print("\n--- Summary ---")
    print(f"  Vs:   {vs.min():.0f} - {vs.max():.0f} m/s")
    print(f"  Vs_f: {vs_filtered_mesh.min():.0f} - {vs_filtered_mesh.max():.0f} m/s")
    print(f"  Vp:   {vp.min():.0f} - {vp.max():.0f} m/s")
    print(f"  Vp_f: {vp_filtered_mesh.min():.0f} - {vp_filtered_mesh.max():.0f} m/s")
    print(f"  dlnVs_filtered: {dlnvs_filtered.min():+.2f} - {dlnvs_filtered.max():+.2f} %")
    print(f"  dlnVp_filtered: {dlnvp_filtered.min():+.2f} - {dlnvp_filtered.max():+.2f} %")
    print("Done.")


if __name__ == "__main__":
    main()
