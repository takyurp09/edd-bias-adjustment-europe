#!/usr/bin/env python

from pathlib import Path
import pickle
import os
import time

import geopandas as gpd
import numpy as np
import xarray as xr

from pc_access import get_items, get_signed_dataset
from edd_calc import calculate_daily_edd

PROJECT_ROOT = Path(__file__).resolve().parent

# -------------------------------------------------------------------
# CONFIG
# -------------------------------------------------------------------

COUNTRY_NAME = os.environ.get("COUNTRY_NAME", "africa")
COUNTRY_TAG  = os.environ.get("COUNTRY_TAG",  "AFR")
ADM_LEVEL    = int(os.environ.get("ADM_LEVEL", "1"))

YEAR_START, YEAR_END = 1994, 2014

source_ids = ["GFDL-ESM4", "MPI-ESM1-2-HR", "UKESM1-0-LL"]
experiment_ids = ["historical"]

DATA_ROOT = Path("/Users/tahmid/Library/CloudStorage/GoogleDrive-tahmid@udel.edu/Other computers/My Laptop/UDel/Taky_research/ci26_biasadj_cmip6/data")

SHAPEFILE      = DATA_ROOT / f"shapefiles/gadm41_{COUNTRY_TAG}_shp/gadm41_{COUNTRY_TAG}_{ADM_LEVEL}.shp"
EDD_OUTPUT_DIR = PROJECT_ROOT / f"output_{COUNTRY_NAME}"
EDD_OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

CALENDAR_DIR    = DATA_ROOT / "ALL_CROPS_netCDF_5min_filled"
CAL_WINDOWS_DIR = PROJECT_ROOT / "output" / COUNTRY_NAME / "calendar_windows"
CAL_WINDOWS_DIR.mkdir(parents=True, exist_ok=True)

# Parquet output dir — same as in ALTERNATIVE_agg_edd.py
PANEL_OUTPUT_DIR = PROJECT_ROOT / "output" / COUNTRY_NAME

# -------------------------------------------------------------------
# BBOX CLAMP
# -------------------------------------------------------------------

CONTINENT_BBOX = {
    "EUR": (-25.0,  27.0,  45.0,  82.0),
    "AFR": (-20.0, -35.0,  52.0,  38.0),
    "ASIA": (25.0, -12.0, 180.0,  82.0),
    "NAM": (-170.0,  7.0,  -52.0, 84.0),
    "SAM": (-82.0, -56.0,  -34.0, 13.0),
    "OCE": (110.0, -50.0,  180.0,  0.0),
}

# -------------------------------------------------------------------
# CALENDARS
# -------------------------------------------------------------------

crop_cals = {
    "wheat":      ["Wheat", "Wheat.Winter"],
    "maize":      ["Maize"],
    "rice":       ["Rice", "Rice.2"],
    "barley":     ["Barley", "Barley.Winter"],
    "potato":     ["Potatoes"],
    "sugar_beet": ["Sugarbeets"],
    "rapeseed":   ["Rapeseed.Winter"],
}

all_calendars = sorted({c for lst in crop_cals.values() for c in lst})

# -------------------------------------------------------------------
# HELPERS
# -------------------------------------------------------------------

def load_country_shapefile():
    gdf = gpd.read_file(SHAPEFILE).to_crs("EPSG:4326")
    minx, miny, maxx, maxy = gdf.total_bounds

    if COUNTRY_TAG in CONTINENT_BBOX:
        clon_min, clat_min, clon_max, clat_max = CONTINENT_BBOX[COUNTRY_TAG]
        minx = max(minx, clon_min)
        miny = max(miny, clat_min)
        maxx = min(maxx, clon_max)
        maxy = min(maxy, clat_max)

    region_bbox = {"lat": slice(miny, maxy), "lon": slice(minx, maxx)}
    return gdf, region_bbox


def _parquet_exists(safe_id: str, year: int) -> bool:
    """Check if the downstream parquet for this model/year already exists."""
    p = PANEL_OUTPUT_DIR / f"panel_edd_{COUNTRY_NAME}_{safe_id}_{year}.parquet"
    return p.exists()


# -------------------------------------------------------------------
# MAIN GENERATOR
# -------------------------------------------------------------------

def generate_edd_files():

    print(f"Reading ADM{ADM_LEVEL} shapefile for {COUNTRY_TAG}...")
    gdf_country, region_bbox = load_country_shapefile()

    items = get_items(source_ids, experiment_ids)

    for item in items:

        if not all(k in item.assets for k in ["tasmax", "tasmin"]):
            continue

        safe_id = str(item.id).replace("/", "_")

        # -------------------------------------------------------
        # KEY FIX: check parquet existence, not .nc existence.
        # .nc files get deleted after processing, so checking them
        # always returns "missing" on a re-run, causing re-download.
        # -------------------------------------------------------
        missing_years = [
            y for y in range(YEAR_START, YEAR_END + 1)
            if not _parquet_exists(safe_id, y)
        ]

        done_years = [
            y for y in range(YEAR_START, YEAR_END + 1)
            if _parquet_exists(safe_id, y)
        ]

        print(f"\n{item.id}")
        print(f"  Already done (parquet exists): {done_years}")
        print(f"  Missing (need to process):     {missing_years}")

        if not missing_years:
            print(f"  ✅ All years done, skipping.")
            continue

        # -------------------------------------------------------
        # Download and compute only missing years
        # -------------------------------------------------------
        print(f"\n  Opening Zarr for {item.id}...")

        for i, year in enumerate(missing_years):

            # Refresh token on every year (i % 1 == 0 is always true,
            # keeping your original intent but now clearly per-year)
            print(f"  🔄 Refreshing token for {item.id} year {year}...")
            try:
                ds_max_full = get_signed_dataset(item, "tasmax", None, region_bbox)
                ds_min_full = get_signed_dataset(item, "tasmin", None, region_bbox)
            except Exception as e:
                print(f"  ERROR refresh {item.id} {year}: {e}")
                continue  # skip this year only, keep going

            print(f"  Processing {item.id} ({year})")

            try:
                ds_tasmax_yr = ds_max_full.sel(
                    time=slice(f"{year}-01-01", f"{year}-12-31")
                )["tasmax"]
                ds_tasmin_yr = ds_min_full.sel(
                    time=slice(f"{year}-01-01", f"{year}-12-31")
                )["tasmin"]
            except Exception as e:
                print(f"  ERROR slice {item.id} {year}: {e}")
                continue

            if ds_tasmax_yr.sizes.get("time", 0) == 0:
                print(f"  WARNING: 0 timesteps for {item.id} {year}, skipping.")
                continue

            try:
                import dask
                ds_tasmax_yr, ds_tasmin_yr = dask.compute(ds_tasmax_yr, ds_tasmin_yr)
            except Exception:
                ds_tasmax_yr = ds_tasmax_yr.load()
                ds_tasmin_yr = ds_tasmin_yr.load()

            ds_tasmax_yr -= 273.15
            ds_tasmin_yr -= 273.15

            edd_ds = calculate_daily_edd(ds_tasmax_yr, ds_tasmin_yr)

            out = EDD_OUTPUT_DIR / f"edd_{COUNTRY_NAME}_{safe_id}_{year}.nc"

            try:
                edd_ds.to_netcdf(out)
                print(f"  Saved {out.name}")
                yield out
            except Exception as e:
                print(f"  ERROR save {item.id} {year}: {e}")

# -------------------------------------------------------------------
# CALENDAR WINDOWS (unchanged)
# -------------------------------------------------------------------

def compute_windows_for_calendar(gdf, cal_name):
    dscal = xr.open_dataset(CALENDAR_DIR / f"{cal_name}.crop.calendar.fill.nc")

    n = len(gdf)

    plant_start = np.full(n, np.nan)
    plant_end   = np.full(n, np.nan)

    for ii, row in gdf.iterrows():
        try:
            pix = dscal.sel(longitude=row.geometry.centroid.x,
                            latitude=row.geometry.centroid.y,
                            method="nearest")

            plant_start[ii] = int(pix["plant.start"])
            plant_end[ii]   = int(pix["plant.end"])
        except:
            continue

    return {"plant": {"start": plant_start, "end": plant_end}}


def precompute_calendar_windows():

    gdf = gpd.read_file(SHAPEFILE).to_crs("EPSG:4326")
    gdf = gdf.reset_index(drop=True)

    for cal in all_calendars:

        out = CAL_WINDOWS_DIR / f"calendar_windows_{COUNTRY_TAG}_{cal}.pkl"

        if out.exists():
            continue

        windows = compute_windows_for_calendar(gdf, cal)

        with open(out, "wb") as f:
            pickle.dump(windows, f)

# -------------------------------------------------------------------
# MAIN
# -------------------------------------------------------------------

if __name__ == "__main__":
    for _ in generate_edd_files():
        pass

    precompute_calendar_windows()