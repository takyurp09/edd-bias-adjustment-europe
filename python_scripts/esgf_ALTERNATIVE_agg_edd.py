#!/usr/bin/env python

import os
import re
import pickle
import time
from pathlib import Path
# from concurrent.futures import ProcessPoolExecutor

import geopandas as gpd
import numpy as np
import pandas as pd
import rasterio
import rioxarray
import xagg as xa
import xarray as xr
import xesmf as xe

xa.set_options(nan_to_zero_regridding=False)

# -------------------------------------------------------------------
# IMPORT COUNTRY CONFIG + EDD/CALENDAR GENERATORS
# -------------------------------------------------------------------

from esgf_ALTERNATIVE_edds_calendar import (
    COUNTRY_NAME,
    COUNTRY_TAG,
    ADM_LEVEL,
    DATA_ROOT,
    crop_cals,                 # calendar names per crop
    generate_edd_files,        # step 1: make EDD .nc
    precompute_calendar_windows,  # step 2: calendar windows pkl
)

# -------------------------------------------------------------------
# CONFIGURATION (derived from imported settings)
# -------------------------------------------------------------------

# Directory with EDD .nc files (from generate_edd_files)
edd_dir = Path(f"output_{COUNTRY_NAME}")

# Output directory for aggregated panel
output_dir = Path("output") / COUNTRY_NAME
output_dir.mkdir(parents=True, exist_ok=True)

# Precomputed calendar windows directory (.pkl files)
calendar_windows_dir = output_dir / "calendar_windows"

# GADM shapefile (ADM_LEVEL)
shapefile = (
    DATA_ROOT
    / f"shapefiles/gadm41_{COUNTRY_TAG}_shp/gadm41_{COUNTRY_TAG}_{ADM_LEVEL}.shp"
)

# Base path for harvested area rasters
harvest_base = DATA_ROOT / "harvested_area_grids"

# EDD variable names in the .nc files
eddlist = ["edd_0", "edd_4", "edd_8", "edd_12", "edd_28", "edd_30", "edd_32"]

# Crop codes (GAEZ / harvested-area codes)
crop_codes = {
    "wheat": 1,
    "maize": 2,
    "rice": 3,
    "barley": 4,
    "potato": 10,
    "sugar_beet": 13,
    "rapeseed": 15
}

country_code = COUNTRY_TAG

# Final merged panel output (all years)
panel_path = output_dir / f"{COUNTRY_NAME}_edd_adm{ADM_LEVEL}_seasonal_panel.parquet"


# -------------------------------------------------------------------
# HELPERS
# -------------------------------------------------------------------

def create_weights(harvestpath: str, gdf: gpd.GeoDataFrame, grid_da: xr.DataArray):
    """
    Create polygon weights for xagg aggregation using a harvested-area raster
    and the grid of `grid_da` (one time slice of the EDD data).
    """
    with rasterio.open(harvestpath) as src:
        data_array = rioxarray.open_rasterio(src)
        weights = data_array.squeeze().rename({"y": "latitude", "x": "longitude"})

        latitudes = weights["latitude"].values
        longitudes = weights["longitude"].values

        delta_lat = np.abs(latitudes[1] - latitudes[0]) / 2
        delta_lon = np.abs(longitudes[1] - longitudes[0]) / 2

        lat_bounds = np.array(
            [[lat - delta_lat, lat + delta_lat] for lat in latitudes]
        )
        lon_bounds = np.array(
            [[lon - delta_lon, lon + delta_lon] for lon in longitudes]
        )

        lat_bounds_da = xr.DataArray(
            lat_bounds,
            dims=["latitude", "bnds"],
            coords={"latitude": latitudes, "bnds": [0, 1]},
        )
        lon_bounds_da = xr.DataArray(
            lon_bounds,
            dims=["longitude", "bnds"],
            coords={"longitude": longitudes, "bnds": [0, 1]},
        )

        weights2 = weights.to_dataset(name="weights")
        weights2["lat_bounds"] = lat_bounds_da
        weights2["lon_bounds"] = lon_bounds_da

        weights2["latitude"].attrs["bounds"] = "lat_bounds"
        weights2["longitude"].attrs["bounds"] = "lon_bounds"

        weights3 = xr.decode_cf(weights2)

    # Regrid harvest map onto the EDD grid
    regridder = xe.Regridder(weights3, grid_da, "bilinear", reuse_weights=False)
    weights4 = regridder(weights3)

    # Pixel overlaps between polygons and regridded harvest weights
    return xa.pixel_overlaps(grid_da, gdf, weights=weights4.weights)


def apply_calendar_windows_vectorized(
    ds: xr.Dataset, windows: dict, varlist: list[str]
) -> xr.Dataset:
    """
    Vectorized calendar-window aggregation.
    ds: Dataset with dims ('poly_idx', 'time') and data_vars in varlist.
    windows: dict loaded from calendar_windows_*.pkl
    varlist: EDD variable names.
    Returns: Dataset with dims ('poly_idx', 'period') and vars in varlist + 'days'.
    """
    # Duplicate time dimension (to handle wrap-around)
    ds2 = xr.concat([ds, ds], dim="time")

    period_order = ["plant", "between", "harvest"]

    t_idx = xr.DataArray(
        np.arange(ds2.sizes["time"]),
        dims=("time",),
        coords={"time": ds2["time"]},
    )

    all_periods = []

    for period in period_order:
        start_arr = windows[period]["start"]  # (n_poly,)
        end_arr   = windows[period]["end"]    # (n_poly,)
        days_arr  = windows[period]["days"]   # (n_poly,)

        start_da = xr.DataArray(start_arr, dims=("poly_idx",))
        end_da   = xr.DataArray(end_arr,   dims=("poly_idx",))
        days_da  = xr.DataArray(days_arr,  dims=("poly_idx",))

        tt, ss = xr.broadcast(t_idx, start_da)
        _, ee  = xr.broadcast(t_idx, end_da)

        mask = (tt >= ss) & (tt < ee)

        masked = ds2[varlist].where(mask, other=0.0)
        summed = masked.sum(dim="time", skipna=True)

        summed = summed.assign_coords(period=period)
        summed["days"] = days_da

        all_periods.append(summed)

    out = xr.concat(all_periods, dim="period")
    out = out.assign_coords(period=("period", period_order))

    return out


# -------------------------------------------------------------------
# WORKER FUNCTION: process a single .nc file
# -------------------------------------------------------------------




#commenting out this cz I wanted to skip the already generated parquet files
"""
def process_nc_file(nc_path_str: str) -> str | None:
   
    Process one EDD .nc file:
      - read shapefile
      - aggregate EDD to ADM polygons with xagg weights
      - apply calendars for wheat/maize/rice (irrigated & rainfed)
      - save one parquet panel per .nc file

    Returns: path to per-file parquet (as str), or None on error.
    nc_path = Path(nc_path_str)
    print(f"Processing {nc_path.name}")
"""


def process_nc_file(nc_path_str: str) -> str | None:
    """
    Process one EDD .nc file:
      - read shapefile
      - aggregate EDD to ADM polygons with xagg weights
      - apply calendars for wheat/maize/rice (irrigated & rainfed)
      - save one parquet panel per .nc file

    Returns: path to per-file parquet (as str), or None on error.
    """
    nc_path = Path(nc_path_str)

    # 🔁 NEW: check if the per-file parquet already exists, and skip if so
    per_file_out = output_dir / f"panel_{nc_path.stem}.parquet"   # NEW
    if per_file_out.exists():                                     # NEW
        print(                                                    # NEW
            f"⚠️  Skipping {nc_path.name} — "                     # NEW
            f"found existing parquet: {per_file_out}"             # NEW
        )                                                         # NEW
        return str(per_file_out)                                  # NEW

    print(f"Processing {nc_path.name}")






    # Read ADM shapefile
    #gdf = gpd.read_file(shapefile, encoding = "latin1")
    gdf = gpd.read_file(shapefile) # ONLY for INDIA
    gdf = gdf.to_crs("EPSG:4326")
    gdf = gdf.reset_index(drop=True)

    adm_codes_col = [c for c in gdf.columns if c.startswith("GID_")][0]
    name_col      = [c for c in gdf.columns if c.startswith("NAME_")][0]

    adm_codes = gdf[adm_codes_col].values
    adm_names = gdf[name_col].values

    # Extract year from filename (last 4 digits before .nc)
    fname = nc_path.name
    year_match = re.search(r"(\d{4})\.nc$", fname)
    year = int(year_match.group(1)) if year_match else None

    # run_id: model+scenario etc.
    if year is not None:
        run_id = fname.replace(f"edd_{COUNTRY_NAME}_", "").replace(f"_{year}.nc", "")
    else:
        run_id = fname

    try:
        ds = xr.open_dataset(nc_path)

        edd_vars_present = [v for v in eddlist if v in ds.data_vars]
        if not edd_vars_present:
            return None

        grid_for_weights = ds[edd_vars_present[0]].isel(time=0, drop=True)

        panel_dfs: list[pd.DataFrame] = []

        # Only wheat, maize, rice
        for crop, code in crop_codes.items():

            for irrigated in [False, True]:
                infix = "IRC" if irrigated else "RFC"
                irrigation_label = "irrigated" if irrigated else "rainfed"

                weightpath = output_dir / f"weights_{crop}_{infix}_{country_code}_{run_id}.wm"

                if weightpath.exists():
                    with open(weightpath, "rb") as fp:
                        weightmap = pickle.load(fp)
                else:
                    harvestpath = (
                        harvest_base
                        / f"ANNUAL_AREA_HARVESTED_{infix}_CROP{code}_HA.ASC"
                    )
                    weightmap = create_weights(str(harvestpath), gdf, grid_for_weights)
                    with open(weightpath, "wb") as fp:
                        pickle.dump(weightmap, fp)

                aggregated = xa.aggregate(ds[edd_vars_present], weightmap)
                regdst2m = aggregated.to_dataset()

                if "feature" in regdst2m.dims and "poly_idx" not in regdst2m.dims:
                    regdst2m = regdst2m.rename({"feature": "poly_idx"})

                if "poly_idx" not in regdst2m.dims:
                    continue

                for cal in crop_cals[crop]:
                    windows_path = (
                        calendar_windows_dir
                        / f"calendar_windows_{country_code}_{cal}.pkl"
                    )
                    if not windows_path.exists():
                        continue

                    with open(windows_path, "rb") as f:
                        cal_windows = pickle.load(f)

                    outdst2m = apply_calendar_windows_vectorized(
                        regdst2m[edd_vars_present], cal_windows, edd_vars_present
                    )

                    outdst2m = outdst2m.reset_coords(drop=True)
                    df = outdst2m.to_dataframe().reset_index()

                    if df.empty:
                        continue

                    df["adm_code"] = adm_codes[df["poly_idx"].values]
                    df["adm_name"] = adm_names[df["poly_idx"].values]

                    df["year"]       = year
                    df["run_id"]     = run_id
                    df["crop"]       = crop
                    df["irrigation"] = irrigation_label
                    df["calendar"]   = cal
                    df["country"]    = country_code

                    panel_dfs.append(df)

        if not panel_dfs:
            return None

        panel = pd.concat(panel_dfs, ignore_index=True)

        per_file_out = output_dir / f"panel_{nc_path.stem}.parquet"
        panel.to_parquet(per_file_out)

        return str(per_file_out)

    except Exception as e:
        print(f"❌ ERROR in process_nc_file({nc_path.name}): {e}")
        raise



# -------------------------------------------------------------------
# MAIN: One year EDD, then agg
# -------------------------------------------------------------------

if __name__ == "__main__":
    start = time.time()

    # Step 1: precompute calendar windows (independent of years/files)
    precompute_calendar_windows()

    # Step 2: generate EDD .nc files one by one, aggregate, and delete each
    per_file_paths: list[Path] = []

    for nc_path in generate_edd_files():
        # nc_path is a pathlib.Path from ALTERNATIVE_edds_calendar.generate_edd_files()
        parquet_path = process_nc_file(str(nc_path))
        if parquet_path is not None:
            per_file_paths.append(Path(parquet_path))

        # Delete the .nc file after processing to save disk space
        try:
            os.remove(nc_path)
            print(f"Deleted {nc_path}")
        except OSError as e:
            print(f"Could not delete {nc_path}: {e}")

    # Step 3: merge all per-file parquet panels into one full panel
    if not per_file_paths:
        print("\nNo per-file panels were created; nothing to merge.")
    else:
        panel_list = [pd.read_parquet(p) for p in per_file_paths]
        full_panel = pd.concat(panel_list, ignore_index=True)
        full_panel.to_parquet(panel_path)
        print(f"\nFull panel shape: {full_panel.shape}")

    print(f"\nTotal runtime: {time.time() - start:.1f} seconds")
