#!/usr/bin/env python3
from __future__ import annotations

import os
import pickle
from pathlib import Path

import geopandas as gpd
import numpy as np
import xarray as xr

from edd_calc import calculate_daily_edd

# -------------------------------------------------------------------
# CONFIGURATION
# -------------------------------------------------------------------

COUNTRY_NAME = os.environ.get("COUNTRY_NAME", "europe")
COUNTRY_TAG  = os.environ.get("COUNTRY_TAG",  "EUR")
ADM_LEVEL    = int(os.environ.get("ADM_LEVEL", "1"))

PROJECT_ROOT = Path(__file__).resolve().parent

DATA_ROOT = Path.home() / (
    "Library/CloudStorage/GoogleDrive-tahmid@udel.edu/"
    "Other computers/My Laptop/UDel/Taky_research/ci26_biasadj_cmip6/data"
)

SHAPEFILE       = DATA_ROOT / f"shapefiles/gadm41_{COUNTRY_TAG}_shp/gadm41_{COUNTRY_TAG}_{ADM_LEVEL}.shp"
EDD_OUTPUT_DIR  = PROJECT_ROOT / f"output_{COUNTRY_NAME}"
CALENDAR_DIR    = DATA_ROOT / "ALL_CROPS_netCDF_5min_filled"
CAL_WINDOWS_DIR = PROJECT_ROOT / "output" / COUNTRY_NAME / "calendar_windows"

EDD_OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
CAL_WINDOWS_DIR.mkdir(parents=True, exist_ok=True)

# -------------------------------------------------------------------
# CROP CALENDARS
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
# STEP 1: ERA5 DAILY tas → EDD FILES
# (data already downloaded; just read, compute EDD, save)
# -------------------------------------------------------------------

def generate_edd_files():
    """
    Generator: for each existing tas_{country}_ERA5_{year}.nc file,
    compute EDD and yield the saved edd_{country}_ERA5_{year}.nc path.
    Skips years where the EDD file already exists.
    """
    tas_files = sorted(EDD_OUTPUT_DIR.glob(f"tas_{COUNTRY_NAME}_ERA5_*.nc"))

    if not tas_files:
        raise FileNotFoundError(
            f"No ERA5 tas files found in {EDD_OUTPUT_DIR}. "
            f"Expected pattern: tas_{COUNTRY_NAME}_ERA5_YYYY.nc"
        )

    for tas_path in tas_files:
        year = tas_path.stem.split("_")[-1]
        edd_path = EDD_OUTPUT_DIR / f"edd_{COUNTRY_NAME}_ERA5_{year}.nc"

        if edd_path.exists():
            print(f"  ⏭  EDD {year} already exists, skipping")
            yield edd_path
            continue

        print(f"  Computing EDD for {year} ...")

        ds = xr.open_dataset(tas_path)

        # ERA5 tas files expected to have tasmax and tasmin already in Celsius
        edd_ds = calculate_daily_edd(
            ds_tasmax = ds["tasmax"],
            ds_tasmin = ds["tasmin"],
        )

        edd_ds.to_netcdf(edd_path)
        print(f"  ✅ saved {edd_path.name}")

        yield edd_path


# -------------------------------------------------------------------
# STEP 2: PRECOMPUTE CALENDAR WINDOWS
# (identical logic to ESGF — uses ALL_CROPS_netCDF_5min_filled)
# -------------------------------------------------------------------

def compute_windows_for_calendar(gdf: gpd.GeoDataFrame, cal_name: str) -> dict:
    dscal_path = CALENDAR_DIR / f"{cal_name}.crop.calendar.fill.nc"
    dscal = xr.open_dataset(str(dscal_path))

    gdf = gdf.reset_index(drop=True)
    n_poly = len(gdf)

    plant_start   = np.full(n_poly, np.nan, dtype=float)
    plant_end     = np.full(n_poly, np.nan, dtype=float)
    between_start = np.full(n_poly, np.nan, dtype=float)
    between_end   = np.full(n_poly, np.nan, dtype=float)
    harvest_start = np.full(n_poly, np.nan, dtype=float)
    harvest_end   = np.full(n_poly, np.nan, dtype=float)
    days_plant    = np.full(n_poly, np.nan, dtype=float)
    days_between  = np.full(n_poly, np.nan, dtype=float)
    days_harvest  = np.full(n_poly, np.nan, dtype=float)

    T      = 365
    offset = T

    for ii, row in gdf.iterrows():
        try:
            pix = dscal.sel(
                longitude = row["centroid"].x,
                latitude  = row["centroid"].y,
                method    = "nearest",
            )

            hs_doy = int(pix["harvest.start"])
            he_doy = int(pix["harvest.end"])
            ps_doy = int(pix["plant.start"])
            pe_doy = int(pix["plant.end"])

            harvest_end_day   = he_doy + offset
            harvest_start_day = hs_doy + offset
            if harvest_start_day > harvest_end_day:
                harvest_start_day -= T

            plant_end_day = pe_doy + offset
            if plant_end_day > harvest_end_day:
                plant_end_day -= T

            plant_start_day = ps_doy + offset
            if plant_start_day > plant_end_day:
                plant_start_day -= T

            between_start_day = plant_end_day + 1
            between_end_day   = harvest_start_day - 1
            if between_end_day < between_start_day:
                between_end_day = between_start_day - 1

            plant_start[ii]   = plant_start_day - 1
            plant_end[ii]     = plant_end_day
            between_start[ii] = between_start_day - 1
            between_end[ii]   = between_end_day
            harvest_start[ii] = harvest_start_day - 1
            harvest_end[ii]   = harvest_end_day

            days_plant[ii]   = plant_end_day   - plant_start_day
            days_between[ii] = between_end_day - between_start_day
            days_harvest[ii] = harvest_end_day - harvest_start_day

        except Exception:
            continue

    return {
        "plant":   {"start": plant_start,   "end": plant_end,   "days": days_plant},
        "between": {"start": between_start, "end": between_end, "days": days_between},
        "harvest": {"start": harvest_start, "end": harvest_end, "days": days_harvest},
    }


def precompute_calendar_windows():
    if not SHAPEFILE.exists():
        raise FileNotFoundError(f"Shapefile not found: {SHAPEFILE}")

    print(f"\nPrecomputing calendar windows for {COUNTRY_TAG}...")
    gdf = gpd.read_file(SHAPEFILE).to_crs("EPSG:4326").reset_index(drop=True)
    gdf["centroid"] = gdf.geometry.centroid

    print(f"Total polygons (ADM{ADM_LEVEL}): {len(gdf)}")
    print(f"Saving calendar windows to: {CAL_WINDOWS_DIR}")

    for cal in all_calendars:
        out_path = CAL_WINDOWS_DIR / f"calendar_windows_{COUNTRY_TAG}_{cal}.pkl"
        if out_path.exists():
            print(f"  {cal}: exists, skipping")
            continue

        windows = compute_windows_for_calendar(gdf, cal)
        with open(out_path, "wb") as f:
            pickle.dump(windows, f)
        print(f"  {cal}: saved")


# -------------------------------------------------------------------
# MAIN
# -------------------------------------------------------------------

if __name__ == "__main__":
    for _path in generate_edd_files():
        pass
    precompute_calendar_windows()
