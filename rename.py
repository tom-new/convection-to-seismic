import xarray as xr
import argparse


def recast_variables(ds: xr.Dataset) -> xr.Dataset:
    """Rename variables in the dataset to match my plotting workflows"""
    coord_mapping = {
        "depth": "r",
        "lon": "lon",
        "lat": "lat",
    }

    var_mapping = {
        "Temperature_K": "T",
        "TemperatureDeviation_K": "dT",
        "dlnVs_percent": "dlnVs_linan_percent",
        "dlnVp_percent": "dlnVp_linan_percent",
        "dlnVs_tofi_percent": "dlnVs_tofi_percent",
        "dlnVp_tofi_percent": "dlnVp_tofi_percent",
    }

    # construct a dictionary of variables to keep, based on the name mapping and what is actually in the dataset
    variables_to_keep = {
        old: new for old, new in var_mapping.items() if old in ds.data_vars
    }

    # subset to just the variables we care about
    ds = ds[variables_to_keep.keys()]

    # rename coordinates and variables
    ds = ds.rename({**coord_mapping, **variables_to_keep})

    # reshape dataset from r, lon, lat to r, lat lon for easier plotting
    ds = ds.transpose("r", "lat", "lon")

    # recast longitude to be in the range [-180, 180) for easier plotting
    ds["lon"] = (ds["lon"] + 180) % 360 - 180

    # reorder longitude to be increasing for easier plotting
    ds = ds.sortby("lon")

    # recast radius to be in m instead of nondimensional model units
    rmax = ds["r"].max().item()
    ds["r"] = ds["r"] * 6371e3 / rmax

    # add depth coordinate for convenience, in km
    ds = ds.assign_coords(
        depth=6371 - ds["r"] / 1e3,
    )

    return ds


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Rename variables in a dataset")
    parser.add_argument("input_file", help="Path to the input .nc file")
    parser.add_argument("output_file", help="Path to the output .nc file")
    args = parser.parse_args()

    ds = xr.open_dataset(args.input_file)
    ds_renamed = recast_variables(ds)
    ds_renamed.to_netcdf(args.output_file)
