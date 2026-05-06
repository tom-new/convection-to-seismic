"""
srts_filter.py — Apply S12RTS, S20RTS, S40RTS tomographic filters to a Vs field.

Input:  single .vtu file (output of convert_to_vs.py), with Vs as point_data.
        Coordinates are non-dimensional Cartesian (RMAX=2.208, D_KM=2891 km).
Output: <stem>_srts_filtered.vtu in the same directory, containing the original
        Vs plus Vs_S40RTS, Vs_S20RTS, Vs_S12RTS.

Pipeline
--------
  mesh pts --IDW--> 181x360 regular grid --> SphericalHarmonicExpansion
                                              expand_batch
                                              reparameterize
                                              S40/S20/S12 filter
                                              evaluate_at_depths
                                              synthesize_batch
  mesh pts <--IDW-- 181x360 regular grid <--

Both IDW weight matrices are precomputed once from the shared horizontal
structure of the extruded mesh. SphericalHarmonicExpansion is built once on
the regular grid; synthesize_batch reuses its Legendre polynomial matrices.

Usage:
    python srts_filter.py <input.vtu>
"""

import sys
import time
from pathlib import Path

import numpy as np
import pyvista as pv
from scipy.spatial import cKDTree

from srts import (
    S12RTS,
    S20RTS,
    S40RTS,
    DepthParameterization,
    SphericalHarmonicExpansion,
)

RMAX = 2.208  # non-dim outer radius
D_KM = 2891.0  # mantle depth in km
LMAX = 40
N_LAT = 181
N_LON = 360
K_IDW = 4
IDW_P = 2
DEPTH_TOL = 1.0  # km — tolerance for grouping mesh points into depth layers


def _tick(label):
    print(f"  {label} ...", flush=True)
    return time.perf_counter()


def _tock(t0):
    print(f"    done ({time.perf_counter() - t0:.1f} s)", flush=True)


def regular_grid(nlat, nlon):
    lons = np.linspace(-180.0, 180.0, nlon, endpoint=False)
    lats = np.linspace(-90.0, 90.0, nlat)
    lon_g, lat_g = np.meshgrid(lons, lats)
    return lon_g.ravel(), lat_g.ravel()


def ll_to_unit_xyz(lon_deg, lat_deg):
    phi = np.radians(lon_deg)
    theta = np.radians(lat_deg)
    return np.stack(
        [
            np.cos(theta) * np.cos(phi),
            np.cos(theta) * np.sin(phi),
            np.sin(theta),
        ],
        axis=1,
    )


def build_idw_weights(src_lon, src_lat, tgt_lon, tgt_lat):
    tree = cKDTree(ll_to_unit_xyz(src_lon, src_lat))
    dist, idx = tree.query(ll_to_unit_xyz(tgt_lon, tgt_lat), k=K_IDW)
    w = 1.0 / np.maximum(dist, 1e-12) ** IDW_P
    w /= w.sum(axis=1, keepdims=True)
    return idx, w


def detect_depth_layers(depth_km, tol=DEPTH_TOL):
    rounded = np.round(depth_km / tol) * tol
    unique = np.sort(np.unique(rounded))
    masks = [rounded == d for d in unique]
    return unique, masks


def depth_boundaries(unique_depths):
    if len(unique_depths) == 1:
        d = unique_depths[0]
        return np.array([max(0.0, d - 1.0), min(2890.0, d + 1.0)])
    mids = 0.5 * (unique_depths[:-1] + unique_depths[1:])
    d0 = max(0.0, unique_depths[0] - (mids[0] - unique_depths[0]))
    d1 = min(2890.0, unique_depths[-1] + (unique_depths[-1] - mids[-1]))
    return np.concatenate([[d0], mids, [d1]])


def idw_apply(weights, idx, source):
    """Vectorised IDW across all layers: source (n_layers, n_src) -> (n_layers, n_tgt)."""
    return (weights * source[:, idx]).sum(axis=2)


def main():
    if len(sys.argv) < 2:
        print(f"Usage: {sys.argv[0]} <input.vtu>")
        sys.exit(1)

    input_path = Path(sys.argv[1])
    output_path = Path(sys.argv[2])

    # ── 1. Load ─────────────────────────────────────────────────────────────
    print(f"Reading {input_path} ...", flush=True)
    mesh = pv.read(input_path)
    xyz = np.asarray(mesh.points)
    vs = np.asarray(mesh.point_data["Vs"], dtype=np.float64)
    n_total = len(vs)
    print(f"  {n_total:,} total points")

    nan_mask = ~np.isfinite(vs)
    if nan_mask.any():
        print(
            f"  WARNING: {nan_mask.sum():,} non-finite Vs values — will replace with layer mean before expansion"
        )

    # ── 2. Spherical coordinates ─────────────────────────────────────────────
    r = np.linalg.norm(xyz, axis=1)
    depth_km = np.clip((RMAX - r) * D_KM, 0.0, 2890.0)
    lat = np.degrees(np.arcsin(np.clip(xyz[:, 2] / r, -1.0, 1.0)))
    lon = np.degrees(np.arctan2(xyz[:, 1], xyz[:, 0]))

    # ── 3. Depth layers ───────────────────────────────────────────────────────
    unique_depths, layer_masks = detect_depth_layers(depth_km)
    n_layers = len(unique_depths)
    n_horiz = layer_masks[0].sum()
    print(
        f"  {n_layers} depth layers, {n_horiz:,} pts/layer  "
        f"({unique_depths[0]:.1f}–{unique_depths[-1]:.1f} km)"
    )

    # Extruded mesh: all layers must share the same horizontal point count.
    assert all(
        m.sum() == n_horiz for m in layer_masks
    ), "Layers have different point counts — non-extruded meshes are not supported."

    # ── 4. Regular intermediate grid ─────────────────────────────────────────
    lon_grid, lat_grid = regular_grid(N_LAT, N_LON)

    mesh_lon_h = lon[layer_masks[0]]
    mesh_lat_h = lat[layer_masks[0]]

    # ── 5. IDW weights — precomputed once ────────────────────────────────────
    t0 = _tick("Forward IDW  (mesh → grid)")
    fwd_idx, fwd_w = build_idw_weights(mesh_lon_h, mesh_lat_h, lon_grid, lat_grid)
    _tock(t0)

    t0 = _tick("Backward IDW (grid → mesh)")
    bwd_idx, bwd_w = build_idw_weights(lon_grid, lat_grid, mesh_lon_h, mesh_lat_h)
    _tock(t0)

    # ── 6. Stack layers, fix NaNs, IDW to grid, subtract per-layer mean ───────
    t0 = _tick("Mapping Vs to regular grid")
    vs_layers = np.stack([vs[m] for m in layer_masks])  # (n_layers, n_horiz)

    if nan_mask.any():
        for i in range(n_layers):
            bad = ~np.isfinite(vs_layers[i])
            if bad.any():
                vs_layers[i, bad] = np.nanmean(vs_layers[i])

    vs_on_grid = idw_apply(fwd_w, fwd_idx, vs_layers)  # (n_layers, n_grid)
    layer_means = vs_on_grid.mean(axis=1)  # (n_layers,)
    vs_on_grid -= layer_means[:, np.newaxis]
    _tock(t0)

    # ── 7. srts pipeline ───────────────────────────────────────────────────────
    t0 = _tick("Building SphericalHarmonicExpansion")
    expander = SphericalHarmonicExpansion(lon_grid, lat_grid, lmax=LMAX)
    _tock(t0)

    t0 = _tick("expand_batch")
    layer_cilms = expander.expand_batch(vs_on_grid)
    _tock(t0)

    t0 = _tick("Reparameterizing")
    projector = DepthParameterization()
    model = projector.reparameterize(list(layer_cilms), depth_boundaries(unique_depths))
    _tock(t0)

    print("  Loading S40RTS, S20RTS, S12RTS ...")
    s40, s20, s12 = S40RTS(), S20RTS(), S12RTS()

    t0 = _tick("Filtering (S40 / S20 / S12)")
    filtered_40 = s40.filter(model)
    sl20 = s20.lmax + 1
    filtered_20 = s20.filter(model[:, :, :sl20, :sl20])
    sl12 = s12.lmax + 1
    filtered_12 = s12.filter(model[:, :, :sl12, :sl12])
    _tock(t0)

    t0 = _tick("Evaluating at layer depths")
    at_depths_40 = DepthParameterization.evaluate_at_depths(filtered_40, unique_depths)
    at_depths_20 = DepthParameterization.evaluate_at_depths(filtered_20, unique_depths)
    at_depths_12 = DepthParameterization.evaluate_at_depths(filtered_12, unique_depths)
    _tock(t0)

    # ── 8. Synthesize back to grid, IDW to mesh, restore mean ────────────────
    t0 = _tick("Synthesizing back to mesh")
    grid_40 = expander.synthesize_batch(at_depths_40)  # (n_layers, n_grid)
    grid_20 = expander.synthesize_batch(at_depths_20)
    grid_12 = expander.synthesize_batch(at_depths_12)

    mesh_40 = idw_apply(bwd_w, bwd_idx, grid_40) + layer_means[:, np.newaxis]
    mesh_20 = idw_apply(bwd_w, bwd_idx, grid_20) + layer_means[:, np.newaxis]
    mesh_12 = idw_apply(bwd_w, bwd_idx, grid_12) + layer_means[:, np.newaxis]
    _tock(t0)

    # Scatter (n_layers, n_horiz) back to flat mesh ordering
    vs_s40 = np.empty(n_total, dtype=np.float64)
    vs_s20 = np.empty(n_total, dtype=np.float64)
    vs_s12 = np.empty(n_total, dtype=np.float64)
    for i, mask in enumerate(layer_masks):
        vs_s40[mask] = mesh_40[i]
        vs_s20[mask] = mesh_20[i]
        vs_s12[mask] = mesh_12[i]

    # ── 9. Write output ───────────────────────────────────────────────────────
    mesh.point_data["Vs_tofi"] = vs_s40
    mesh.point_data["Vs_S20RTS_tofi"] = vs_s20
    mesh.point_data["Vs_S12RTS_tofi"] = vs_s12
    mesh.save(output_path)

    print(f"\nDone. Written to {output_path}")
    print(f"  Vs_tofi: {vs_s40.min():.0f} – {vs_s40.max():.0f} m/s")
    print(f"  Vs_S20RTS_tofi: {vs_s20.min():.0f} – {vs_s20.max():.0f} m/s")
    print(f"  Vs_S12RTS_tofi: {vs_s12.min():.0f} – {vs_s12.max():.0f} m/s")


if __name__ == "__main__":
    main()
