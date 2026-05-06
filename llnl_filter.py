"""
tofi_filter.py — Tomographic filtering of Vs and Vp using the
LLNL-G3D-JPS resolution matrix with IDW interpolation.

Usage:
    python tofi_filter.py <input.vtu> <output.vtu>

The input is the output of convert_to_vs.py: a VTU (or PVTU) with Vs and Vp
as point-data fields (m/s) and non-dimensional Cartesian coordinates
(RMAX=2.208, D_KM=2891).  The script adds Vs_filtered and Vp_filtered
fields and writes a new VTU.

Forward projection (mesh → LLNL grid) uses the layer-by-layer IDW in
llnltofi.interpolation.project_onto_grid.  Back-projection (LLNL → mesh)
uses a 3-D inverse-distance-weighted average over the k nearest LLNL nodes.
"""

import sys
from pathlib import Path

import numpy as np
import pyvista as pv
from scipy.spatial import cKDTree

import llnltofi
from llnltofi.interpolation import project_onto_grid
from llnltofi._constants import (
    N_LAYERS_UM_TZ,
    N_LAYERS_LM,
    N_POINTS_UM_TZ,
    N_POINTS_LM,
    R_EARTH_KM,
)

RMAX = 2.208
D_KM = 2891.0
IDW_K_BACK = 8  # LLNL neighbours used for back-projection IDW


# ── coordinate helpers ────────────────────────────────────────────────────────


def nondim_to_spherical(coords):
    """Non-dimensional Cartesian → (gc_lat_deg, lon_deg, radius_km).

    This is the format expected by project_onto_grid.
    """
    r = np.linalg.norm(coords, axis=1)
    depth_km = (RMAX - r) * D_KM
    radius_km = R_EARTH_KM - depth_km
    unit = coords / r[:, np.newaxis]
    lat_deg = np.degrees(np.arcsin(np.clip(unit[:, 2], -1.0, 1.0)))
    lon_deg = np.degrees(np.arctan2(unit[:, 1], unit[:, 0]))
    return np.column_stack([lat_deg, lon_deg, radius_km])


def nondim_to_physical_xyz(coords):
    """Non-dimensional Cartesian → physical Cartesian in metres."""
    r = np.linalg.norm(coords, axis=1)
    depth_km = (RMAX - r) * D_KM
    radius_m = (R_EARTH_KM - depth_km) * 1e3
    unit = coords / r[:, np.newaxis]
    return unit * radius_m[:, np.newaxis]


# ── filtering helpers ─────────────────────────────────────────────────────────


def layer_mean_1d(slowness):
    """Depth-dependent 1D reference from layer-wise mean slowness."""
    n_um_tz = N_LAYERS_UM_TZ * N_POINTS_UM_TZ
    s_um = slowness[:n_um_tz].reshape(N_LAYERS_UM_TZ, N_POINTS_UM_TZ)
    s_lm = slowness[n_um_tz:].reshape(N_LAYERS_LM, N_POINTS_LM)
    return np.concatenate(
        [
            np.repeat(s_um.mean(axis=1), N_POINTS_UM_TZ),
            np.repeat(s_lm.mean(axis=1), N_POINTS_LM),
        ]
    )


def apply_filter(velocity_on_llnl, model):
    """Filter one velocity field through the resolution matrix.

    Converts to slowness, removes the layer-mean 1D reference, applies R,
    then recovers velocity.  No amplitude scaling is needed for either Vp
    or Vs because R operates on the actual slowness anomalies regardless of
    wavetype (see Simmons et al. 2019, Sec. 3).
    """
    s = 1.0 / velocity_on_llnl
    s_1d = layer_mean_1d(s)
    model.values = s - s_1d
    ds_filtered = model.apply()
    return 1.0 / (s_1d + ds_filtered)


# ── back-projection ───────────────────────────────────────────────────────────


def back_project_idw(llnl_values, llnl_xyz, mesh_xyz, k=IDW_K_BACK):
    """IDW back-projection from LLNL grid → arbitrary mesh points.

    Finds the k nearest LLNL grid nodes for every mesh point in 3-D
    physical space (metres) and returns the inverse-distance-weighted mean.
    """
    tree = cKDTree(llnl_xyz)
    dists, idx = tree.query(mesh_xyz, k=k, workers=1)
    dists = np.where(dists == 0.0, 1e-10, dists)  # guard exact coincidences
    weights = 1.0 / dists
    weights /= weights.sum(axis=1, keepdims=True)
    return np.einsum("ij,ij->i", weights, llnl_values[idx])


# ── main ──────────────────────────────────────────────────────────────────────


def main():
    if len(sys.argv) < 3:
        print(f"Usage: {sys.argv[0]} <input.vtu> <output.vtu>")
        sys.exit(1)

    input_path = sys.argv[1]
    output_path = sys.argv[2]

    # ── 1. Load LLNL model ─────────────────────────────────────────────────
    print("Loading LLNL resolution model...")
    model = llnltofi.ResolutionModel()
    llnl_xyz = model.coordinates_in_xyz  # (N_MODEL, 3) metres
    _ = model.R  # ensure R.npz is present
    print(f"  {model.n_model:,} grid points, {model.n_layers} layers")

    # ── 2. Load simulation mesh ────────────────────────────────────────────
    print(f"\nReading {input_path} ...")
    mesh = pv.read(input_path)
    coords = np.asarray(mesh.points)
    vs = np.asarray(mesh.point_data["Vs"])
    vp = np.asarray(mesh.point_data["Vp"])
    print(f"  {len(vs):,} mesh points")

    # ── 3. Convert mesh coordinates ────────────────────────────────────────
    print("Converting coordinates...")
    sph_coords = nondim_to_spherical(coords)  # (N, 3): gc_lat_deg, lon_deg, radius_km
    mesh_xyz = nondim_to_physical_xyz(coords)  # (N, 3): metres

    # ── 4. Forward IDW: mesh → LLNL grid ──────────────────────────────────
    print("\nForward IDW: Vs mesh → LLNL grid...")
    vs_on_llnl = project_onto_grid(sph_coords, vs, model)
    print("Forward IDW: Vp mesh → LLNL grid...")
    vp_on_llnl = project_onto_grid(sph_coords, vp, model)

    # ── 5. Apply resolution filter ─────────────────────────────────────────
    print("\nFiltering Vs through resolution matrix...")
    vs_filtered_llnl = apply_filter(vs_on_llnl, model)
    print("Filtering Vp through resolution matrix...")
    vp_filtered_llnl = apply_filter(vp_on_llnl, model)

    # ── 6. Backward IDW: LLNL grid → mesh ─────────────────────────────────
    print(f"\nBack-projection IDW (k={IDW_K_BACK}): Vs LLNL → mesh...")
    vs_filtered_mesh = back_project_idw(vs_filtered_llnl, llnl_xyz, mesh_xyz)
    print(f"Back-projection IDW (k={IDW_K_BACK}): Vp LLNL → mesh...")
    vp_filtered_mesh = back_project_idw(vp_filtered_llnl, llnl_xyz, mesh_xyz)

    # ── 7. Write output ────────────────────────────────────────────────────
    mesh.point_data["Vs_tofi"] = vs_filtered_mesh
    mesh.point_data["Vp_tofi"] = vp_filtered_mesh

    Path(output_path).parent.mkdir(parents=True, exist_ok=True)
    print(f"\nWriting {output_path} ...")
    mesh.save(output_path)

    print("\n--- Summary ---")
    print(f"  Vs:   {vs.min():.0f} – {vs.max():.0f} m/s")
    print(f"  Vs_tofi: {vs_filtered_mesh.min():.0f} – {vs_filtered_mesh.max():.0f} m/s")
    print(f"  Vp:   {vp.min():.0f} – {vp.max():.0f} m/s")
    print(f"  Vp_tofi: {vp_filtered_mesh.min():.0f} – {vp_filtered_mesh.max():.0f} m/s")
    print("Done.")


if __name__ == "__main__":
    main()
