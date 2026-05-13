"""
Convert non-dimensional mantle convection output (PVTU) to seismic velocities.

Loads the full dataset at once, dimensionalises temperature, then computes Vs
and Vp using the SLB_24 pyroliteCFMASNaCr thermodynamic model with Cammarano Q3
anelastic correction. The thermodynamic table is regularised against the
spherically-averaged temperature profile extracted from the mesh itself, which
removes spurious phase-transition jumps in the Vs(T) relationship. Output
contains only Temperature_K, Vs, and Vp.

Usage:
    python convert_to_vs.py <input.pvtu> <output.vtu>
"""

import sys
import numpy as np
import pyvista as pv
import gdrift

from _layer_mean import dln_percent_by_layer

RMAX = 2.208  # non-dim outer radius (surface)
D_KM = 2891.0  # mantle depth in km
TS = 300.0  # surface temperature (K)
DELTA_T = 3700.0  # temperature drop across mantle (K)
TEMP_FIELD = "FullTemperature_CG"
D_TEMP_FIELD = "Temperature_Deviation_CG"
Q_PROFILE = "Q6"
N_BINS = 200  # probe depths for radial temperature average


def build_average_temperature_profile(depth_m, t_kelvin, n_bins=N_BINS):
    """Compute a spherically-averaged temperature profile from mesh data.

    Divides the mantle into n_bins evenly-spaced depth intervals and averages
    the temperature of all mesh points that fall within each interval. Bins
    with no points are silently dropped, so the returned arrays may be shorter
    than n_bins if the mesh resolution is coarser than expected.

    Parameters
    ----------
    depth_m : ndarray
        Depth of every mesh node in metres (surface = 0, CMB = D_KM*1e3).
    t_kelvin : ndarray
        Temperature of every mesh node in Kelvin.
    n_bins : int
        Number of evenly-spaced probe depths to attempt.

    Returns
    -------
    gdrift.SplineProfile
        Radial temperature profile suitable for use as the regularisation
        anchor in regularise_thermodynamic_table.
    """
    bin_edges = np.linspace(0, D_KM * 1e3, n_bins + 1)
    bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])

    valid_depths, avg_temps = [], []
    for lo, hi, center in zip(bin_edges[:-1], bin_edges[1:], bin_centers):
        mask = (depth_m >= lo) & (depth_m < hi)
        if mask.any():
            valid_depths.append(center)
            avg_temps.append(t_kelvin[mask].mean())

    n_valid = len(valid_depths)
    print(f"  Radial average: {n_valid}/{n_bins} bins populated")

    # Prepend depth=0 with the known surface temperature so the profile
    # covers the full table range and regularise_thermodynamic_table
    # does not produce NaN for the shallowest mesh layer.
    valid_depths.insert(0, 0.0)
    avg_temps.insert(0, TS)

    return gdrift.SplineProfile(
        depth=np.array(valid_depths),
        value=np.array(avg_temps),
        extrapolate=True,
    )


def main():
    if len(sys.argv) < 3:
        print(f"Usage: {sys.argv[0]} <input.pvtu> <output.vtu>")
        sys.exit(1)

    input_pvtu = sys.argv[1]
    output_vtu = sys.argv[2]

    print(f"Reading {input_pvtu} ...")
    mesh = pv.read(input_pvtu)

    coords = np.asarray(mesh.points)
    depth_m = (RMAX - np.linalg.norm(coords, axis=1)) * D_KM * 1e3
    t_kelvin = TS + DELTA_T * np.asarray(mesh.point_data[TEMP_FIELD])
    dt_kelvin = mesh.point_data[D_TEMP_FIELD] * DELTA_T

    print("Building radial temperature profile ...")
    temperature_profile = build_average_temperature_profile(depth_m, t_kelvin)

    print("Loading SLB_24 pyroliteCFMASNaCr ...")
    slb = gdrift.ThermodynamicModel("SLB_24", "pyroliteCFMASNaCr")

    print("Regularising thermodynamic table ...")
    regular_slb = gdrift.regularise_thermodynamic_table(
        slb,
        temperature_profile,
        regular_range={
            "v_s": (-1.5, 0.0),
            "v_p": (-np.inf, 0.0),
            "rho": (-np.inf, 0.0),
        },
    )

    print(f"Applying Cammarano {Q_PROFILE} anelastic correction ...")
    anelastic = gdrift.CammaranoAnelasticityModel.from_q_profile(Q_PROFILE)
    corrected = gdrift.apply_anelastic_correction(regular_slb, anelastic)

    depth_min = slb.get_depths().min()
    depth_max = slb.get_depths().max()
    temp_min = slb.get_temperatures().min()
    temp_max = slb.get_temperatures().max()

    depth_c = np.clip(depth_m, depth_min, depth_max)
    temp_c = np.clip(t_kelvin, temp_min, temp_max)

    print("Computing Vs and Vp ...")
    vs = corrected.temperature_to_vs(temp_c, depth_c)
    vp = corrected.temperature_to_vp(temp_c, depth_c)

    print("Computing dlnVs and dlnVp ...")
    depth_km = depth_m / 1e3
    dlnvs = dln_percent_by_layer(vs, depth_km)
    dlnvp = dln_percent_by_layer(vp, depth_km)

    out_mesh = pv.UnstructuredGrid(mesh.cells, mesh.celltypes, mesh.points)
    out_mesh.point_data["T"] = t_kelvin
    out_mesh.point_data["dT"] = dt_kelvin
    out_mesh.point_data["Vs"] = vs
    out_mesh.point_data["Vp"] = vp
    out_mesh.point_data["dlnVs"] = dlnvs
    out_mesh.point_data["dlnVp"] = dlnvp

    print(f"Writing {output_vtu} ...")
    out_mesh.save(output_vtu)

    print("\n--- Summary ---")
    print(f"Points:      {len(t_kelvin):,}")
    print(f"Depth range: {depth_m.min()/1e3:.1f} - {depth_m.max()/1e3:.1f} km")
    print(f"T range:     {t_kelvin.min():.0f} - {t_kelvin.max():.0f} K")
    print(f"dT range:    {dt_kelvin.min():.0f} - {dt_kelvin.max():.0f} K")
    print(f"Vs range:    {np.nanmin(vs):.0f} - {np.nanmax(vs):.0f} m/s")
    print(f"Vp range:    {np.nanmin(vp):.0f} - {np.nanmax(vp):.0f} m/s")
    print(f"dlnVs range: {np.nanmin(dlnvs):+.2f} - {np.nanmax(dlnvs):+.2f} %")
    print(f"dlnVp range: {np.nanmin(dlnvp):+.2f} - {np.nanmax(dlnvp):+.2f} %")
    print("Done.")


if __name__ == "__main__":
    main()
