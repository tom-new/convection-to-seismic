"""
Convert filtered seismic velocities (Vs, Vp) to deviation from layer mean
(dlnVs, dlnVp).

Loads the full netCDF dataset at once, calculates area weighted radial averages
of Vs and Vp, then computes dlnVs and dlnVp as (V - V_avg) / V_avg. Output
contains Temperature_K, Vs, Vp, dlnVs, and dlnVp.

Usage:
    python convert_to_dlnV.py <input.nc> <output.nc>
"""

import sys
import numpy as np
import xarray as xr

field_mapping = {
    "Vs": "dlnVs",
    "Vp": "dlnVp",
    "Vs_filtered": "dlnVs_tofi",
    "Vp_filtered": "dlnVp_tofi",
}


def build_average_profile(da):
    """Compute a spherically-averaged profile from a DataArray.

    Divides the mantle into n_bins evenly-spaced depth intervals and averages
    the temperature of all mesh points that fall within each interval. Bins
    with no points are silently dropped, so the returned arrays may be shorter
    than n_bins if the mesh resolution is coarser than expected.

    Parameters
    ----------
    da : xarray.DataArray
        DataArray containing the values to be averaged. Must have coordinates (r, lat, lon).

    Returns
    -------
    xarray.DataArray
        Horizontally averaged DataArray with the same coordinates and shape as the input (i.e. still 3D but with the values replaced by the average at each depth).
    """

    lat_weights = np.cos(
        np.radians(da.lat)
    )  # convert to radians then compute cos(lat) to get spherical surface area weights
    lat_weights = lat_weights.clip(min=0)  # avoid tiny negative weights near the poles
    mean_by_r = da.mean(
        dim=("lon", "lat")
    )  # compute the weighted mean along the lat and lon dimensions, leaving only the r dimension

    return mean_by_r


def process_variable(da: xr.DataArray) -> xr.DataArray:
    """Compute the dlnV version of a velocity DataArray."""
    avg_profile = build_average_profile(da)
    dlnV = (da - avg_profile) / avg_profile
    return dlnV


def main():
    if len(sys.argv) < 2 or len(sys.argv) > 3:
        print(f"Usage: {sys.argv[0]} <input.nc> [<output.nc>]")
        print(
            "  <input.nc>   - path to the input netCDF file containing Temperature_K, Vs, and Vp"
        )
        print(
            "  <output.nc>  - path to the output netCDF file to write (default: overwrite input)"
        )
        sys.exit(1)

    input_path = sys.argv[1]
    output_path = sys.argv[2] if len(sys.argv) == 3 else input_path
    print(f"Reading {input_path} ...")
    ds = xr.open_dataset(input_path)

    present_vars = {
        old_name: new_name
        for old_name, new_name in field_mapping.items()
        if old_name in ds.data_vars
    }

    # set zero values to NaN so they are ignored in future calculations
    for var in ds.data_vars:
        ds[var] = ds[var].where(ds[var] != 0, np.nan)

    print("Computing dlnVs and dlnVp ...")
    for old_name, new_name in present_vars.items():
        print(f"... Processing {old_name} → {new_name} ...")
        ds[new_name] = process_variable(ds[old_name])
        ds[new_name + "_percent"] = ds[new_name] * 100

    print(f"Writing {output_path} ...")
    ds.to_netcdf(output_path)

    print("\n--- Summary ---")
    if "Vs" in ds.data_vars:
        print(
            f"Vs range:       {np.nanmin(ds['Vs']):.0f} - {np.nanmax(ds['Vs']):.0f} m/s"
        )
    if "Vp" in ds.data_vars:
        print(
            f"Vp range:       {np.nanmin(ds['Vp']):.0f} - {np.nanmax(ds['Vp']):.0f} m/s"
        )
    if "dlnVs_percent" in ds.data_vars:
        print(
            f"dlnVs range:    {np.nanmin(ds['dlnVs_percent']):.1f} - {np.nanmax(ds['dlnVs_percent']):.2f}%"
        )
    if "dlnVp_percent" in ds.data_vars:
        print(
            f"dlnVp range:    {np.nanmin(ds['dlnVp_percent']):.1f} - {np.nanmax(ds['dlnVp_percent']):.2f}%"
        )
    print("Done.")


if __name__ == "__main__":
    main()
