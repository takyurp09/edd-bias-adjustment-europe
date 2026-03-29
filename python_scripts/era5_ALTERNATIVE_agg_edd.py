#!/usr/bin/env python

import os
import re
import pickle
import time
from pathlib import Path

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

from ALTERNATIVE_edds_calendar import (
    COUNTRY_NAME,
    COUNTRY_TAG,
    ADM_LEVEL,
    DATA_ROOT,
    EDD_OUTPUT_DIR,
    CAL_WINDOWS_DIR,
    crop_cals,
    generate_edd_files,
    precompute_calendar_windows,
)

# -------------------------------------------------------------------
# CONFIGURATION
# -------------------------------------------------------------------

# Output directory for aggregated panel
output_dir = Path("output") / COUNTRY_NAME
output_dir.mkdir(parents=True, exist_ok=True)

# GADM shapefile
shapefile = (
    DATA_ROOT
    / f"shapefiles/gadm41_{COUNTRY_TAG}_shp/gadm41_{COUNTRY_TAG}_{ADM_LEVEL}.shp"
)

# Harvested area rasters
harvest_base = DATA_ROOT / "harvested_area_grids"

# Expected EDD variables
eddlist = ["edd_0", "edd_4", "edd_8", "edd_12", "edd_28", "edd_30", "edd_32"]

# Crop codes (GAEZ / harvested-area codes)
crop_codes = {
    "wheat":      1,
    "maize":      2,
    "rice":       3,
    "barley":     4,
    "potato":     10,
    "sugar_beet": 13,
    "rapeseed":   15,
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
    and the grid of grid_da (one time slice of the EDD data).
    ERA5 grid is fixed across years — weights are reused if .wm already exists.
    """
    with rasterio.open(harvestpath) as src:
        data_array = rioxarray.open_rasterio(src)
        weights = data_array.squeeze().rename({"y": "latitude", "x": "longitude"})

        latitudes  = weights["latitude"].values
        longitudes = weights["longitude"].values

        delta_lat = np.abs(latitudes[1]  - latitudes[0]) / 2
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

        weights2["latitude"].attrs["bounds"]  = "lat_bounds"
        weights2["longitude"].attrs["bounds"] = "lon_bounds"

        weights3 = xr.decode_cf(weights2)

    regridder = xe.Regridder(weights3, grid_da, "bilinear", reuse_weights=False)
    weights4  = regridder(weights3)

    return xa.pixel_overlaps(grid_da, gdf, weights=weights4.weights)


def apply_calendar_windows_vectorized(
    ds: xr.Dataset, windows: dict, varlist: list
) -> xr.Dataset:
    """
    Vectorized calendar-window aggregation.
    Identical to ESGF version.
    """
    ds2 = xr.concat([ds, ds], dim="time")

    period_order = ["plant", "between", "harvest"]

    t_idx = xr.DataArray(
        np.arange(ds2.sizes["time"]),
        dims=("time",),
        coords={"time": ds2["time"]},
    )

    all_periods = []

    for period in period_order:
        start_da = xr.DataArray(windows[period]["start"], dims=("poly_idx",))
        end_da   = xr.DataArray(windows[period]["end"],   dims=("poly_idx",))
        days_da  = xr.DataArray(windows[period]["days"],  dims=("poly_idx",))

        tt, ss = xr.broadcast(t_idx, start_da)
        _, ee  = xr.broadcast(t_idx, end_da)

        mask   = (tt >= ss) & (tt < ee)
        summed = ds2[varlist].where(mask, other=0.0).sum(dim="time", skipna=True)

        summed = summed.assign_coords(period=period)
        summed["days"] = days_da

        all_periods.append(summed)

    out = xr.concat(all_periods, dim="period")
    out = out.assign_coords(period=("period", period_order))

    return out


# -------------------------------------------------------------------
# PROCESS ONE EDD FILE
# -------------------------------------------------------------------

def process_nc_file(nc_path_str: str):
    """
    Process one EDD .nc file:
      - aggregate EDD to ADM polygons with xagg weights
      - apply calendar windows for all crops (irrigated & rainfed)
      - save one parquet panel per .nc file

    Key ERA5 difference vs ESGF:
      - .wm weight files have NO run_id suffix (ERA5 grid is fixed across years)
      - no model loop, no checkpoint — just one ERA5 source per year
    """
    nc_path = Path(nc_path_str)

    per_file_out = output_dir / f"panel_{nc_path.stem}.parquet"
    if per_file_out.exists():
        print(f"  ⏭  Skipping {nc_path.name} — existing parquet: {per_file_out}")
        return str(per_file_out)

    print(f"Processing {nc_path.name}")

    gdf = gpd.read_file(shapefile).to_crs("EPSG:4326").reset_index(drop=True)
    adm_code_col = [c for c in gdf.columns if c.startswith("GID_")][0]
    adm_name_col = [c for c in gdf.columns if c.startswith("NAME_")][0]

    adm_codes = gdf[adm_code_col].values
    adm_names = gdf[adm_name_col].values

    year_match = re.search(r"(\d{4})\.nc$", nc_path.name)
    year   = int(year_match.group(1)) if year_match else None
    run_id = f"ERA5_{year}"

    try:
        ds = xr.open_dataset(nc_path)
        vars_present = [v for v in eddlist if v in ds.data_vars]

        if not vars_present:
            return None

        grid_for_weights = ds[vars_present[0]].isel(time=0, drop=True)

        panel_dfs = []

        for crop, code in crop_codes.items():
            for irrigated in [False, True]:
                infix            = "IRC" if irrigated else "RFC"
                irrigation_label = "irrigated" if irrigated else "rainfed"

                # ERA5 grid is fixed — share weights across years (no run_id in name)
                weightpath = output_dir / f"weights_{crop}_{infix}_{country_code}.wm"

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

                aggregated = xa.aggregate(ds[vars_present], weightmap).to_dataset()

                if "feature" in aggregated.dims and "poly_idx" not in aggregated.dims:
                    aggregated = aggregated.rename({"feature": "poly_idx"})

                if "poly_idx" not in aggregated.dims:
                    continue

                for cal in crop_cals[crop]:
                    windows_path = (
                        CAL_WINDOWS_DIR
                        / f"calendar_windows_{country_code}_{cal}.pkl"
                    )
                    if not windows_path.exists():
                        continue

                    with open(windows_path, "rb") as f:
                        cal_windows = pickle.load(f)

                    out = apply_calendar_windows_vectorized(
                        aggregated[vars_present], cal_windows, vars_present
                    )

                    out = out.reset_coords(drop=True)
                    df  = out.to_dataframe().reset_index()

                    if df.empty:
                        continue

                    df["adm_code"]  = adm_codes[df["poly_idx"].values]
                    df["adm_name"]  = adm_names[df["poly_idx"].values]
                    df["year"]      = year
                    df["run_id"]    = run_id
                    df["crop"]      = crop
                    df["irrigation"] = irrigation_label
                    df["calendar"]  = cal
                    df["country"]   = country_code

                    panel_dfs.append(df)

        if not panel_dfs:
            return None

        panel = pd.concat(panel_dfs, ignore_index=True)
        panel.to_parquet(per_file_out)

        return str(per_file_out)

    except Exception as e:
        print(f"❌ ERROR in process_nc_file({nc_path.name}): {e}")
        raise


# -------------------------------------------------------------------
# MAIN
# -------------------------------------------------------------------

if __name__ == "__main__":
    start = time.time()

    # Step 1: precompute calendar windows (once, shared across all years)
    precompute_calendar_windows()

    # Step 2: generate EDD .nc files from existing ERA5 tas files, then aggregate
    # NOTE: unlike ESGF, we do NOT delete .nc files after processing
    # (ERA5 data is pre-downloaded and kept on disk)
    per_file_paths = []

    for nc_path in generate_edd_files():
        parquet_path = process_nc_file(str(nc_path))
        if parquet_path is not None:
            per_file_paths.append(Path(parquet_path))

    # Step 3: merge all per-year parquet panels into one full panel
    if not per_file_paths:
        print("\nNo per-file panels were created; nothing to merge.")
    else:
        full_panel = pd.concat(
            [pd.read_parquet(p) for p in per_file_paths], ignore_index=True
        )
        full_panel.to_parquet(panel_path)
        print(f"\nFull panel shape: {full_panel.shape}")
        print(f"WROTE: {panel_path}")

    print(f"\nTotal runtime: {time.time() - start:.1f}s")
