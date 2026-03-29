#!/usr/bin/env python

from pathlib import Path
import pickle
import os

import geopandas as gpd
import numpy as np
import xarray as xr

from pangeo_access import get_items, get_signed_dataset


from esgf_edd_calc import calculate_daily_edd  # daily EDD calculator

PROJECT_ROOT = Path(__file__).resolve().parent


# -------------------------------------------------------------------
# CONFIGURATION
# -------------------------------------------------------------------

# Which country? (read from env vars, with defaults)
COUNTRY_NAME = os.environ.get("COUNTRY_NAME", "europe")
COUNTRY_TAG  = os.environ.get("COUNTRY_TAG", "EUR")
ADM_LEVEL    = int(os.environ.get("ADM_LEVEL", "1"))

def lon360_to_lon180(da):
    if "lon" in da.coords and float(da["lon"].max()) > 180:
        da = da.assign_coords(lon=((da["lon"] + 180) % 360) - 180).sortby("lon")
    return da

# Year range
YEAR_START, YEAR_END = 1994, 2014

# 5 Europe-focused models not in GDPCIR
#source_ids = [
 #   "IPSL-CM6A-LR",    # French — excellent European temperature/precip
  #  "CNRM-CM6-1",      # French/Météo-France — best for European climate
   # "CNRM-CM6-1-HR",   # Same but HIGH resolution — great for Europe
    #"EC-Earth3",       # European consortium — built for European climate
    #"MRI-ESM2-0",      # Top performer over Europe in CMIP6 evaluations
#]


source_ids = ["GFDL-ESM4", "MPI-ESM1-2-HR"]

experiment_ids = ["historical"]

# Root data directory
DATA_ROOT = Path.home() / "Library/CloudStorage/GoogleDrive-tahmid@udel.edu/Other computers/My Laptop/UDel/Taky_research/ci26_biasadj_cmip6/data"

# Paths
SHAPEFILE      = DATA_ROOT / f"shapefiles/gadm41_{COUNTRY_TAG}_shp/gadm41_{COUNTRY_TAG}_{ADM_LEVEL}.shp"
EDD_OUTPUT_DIR = PROJECT_ROOT / f"output_{COUNTRY_NAME}"
CALENDAR_DIR   = DATA_ROOT / "ALL_CROPS_netCDF_5min_filled"
CAL_WINDOWS_DIR = PROJECT_ROOT / "output" / COUNTRY_NAME / "calendar_windows"

EDD_OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
CAL_WINDOWS_DIR.mkdir(parents=True, exist_ok=True)

# Checkpoint file — tracks completed (model, year) pairs
CHECKPOINT_FILE = PROJECT_ROOT / f"output_{COUNTRY_NAME}" / "esgf_checkpoint.txt"

# Crop calendars
crop_cals = {
    "wheat":     ["Wheat", "Wheat.Winter"],
    "maize":     ["Maize"],
    "rice":      ["Rice", "Rice.2"],
    "barley":    ["Barley", "Barley.Winter"],
    "potato":    ["Potatoes"],
    "sugar_beet":["Sugarbeets"],
    "rapeseed":  ["Rapeseed.Winter"],
}
all_calendars = sorted({c for lst in crop_cals.values() for c in lst})


# -------------------------------------------------------------------
# CHECKPOINT HELPERS
# -------------------------------------------------------------------

def load_checkpoint() -> set:
    """Load set of completed 'source_id|year' keys."""
    if not CHECKPOINT_FILE.exists():
        return set()
    with open(CHECKPOINT_FILE) as f:
        return set(line.strip() for line in f if line.strip())

def save_checkpoint(source_id: str, year: int):
    """Append a completed key to the checkpoint file."""
    with open(CHECKPOINT_FILE, "a") as f:
        f.write(f"{source_id}|{year}\n")

def checkpoint_key(source_id: str, year: int) -> str:
    return f"{source_id}|{year}"


# -------------------------------------------------------------------
# COMMON: read shapefile and get bbox
# -------------------------------------------------------------------

def load_country_shapefile():
    gdf = gpd.read_file(SHAPEFILE)
    gdf = gdf.to_crs("EPSG:4326")
    minx, miny, maxx, maxy = gdf.total_bounds
    region_bbox = {"lat": slice(miny, maxy), "lon": slice(minx, maxx)}
    return gdf, region_bbox


# -------------------------------------------------------------------
# PART 1: Generate EDD .nc files (per year × model)
# -------------------------------------------------------------------

def generate_edd_files():
    """
    Generator: yields each saved EDD NetCDF file (Path) as it is created.
    Skips already-completed (model, year) pairs via checkpoint file.
    Skips models that fail gracefully so the run continues.
    """
    print(f"Reading ADM{ADM_LEVEL} shapefile for {COUNTRY_TAG}...")
    gdf_country, region_bbox = load_country_shapefile()
    print(f"Total polygons (ADM{ADM_LEVEL}): {len(gdf_country)}")

    completed = load_checkpoint()
    print(f"Checkpoint: {len(completed)} (model, year) pairs already done")

    items = get_items(source_ids, experiment_ids)
    print(f"Found {len(items)} items from ESGF")

    for year in range(YEAR_START, YEAR_END + 1):
        os.environ["CMIP6_YEAR"] = str(year)
        print(f"\n=== YEAR {year} ===")

        for item in items:
            if not all(k in item.assets for k in ["tasmax", "tasmin"]):
                continue

            ck = checkpoint_key(item.source_id, year)

            # --- Resume: skip if already done ---
            if ck in completed:
                print(f"  ⏭  {item.id} ({year}) — already done, skipping")
                # Still yield the existing file if it exists
                safe_id = str(item.id).replace("/", "_")
                edd_outfile = EDD_OUTPUT_DIR / f"edd_{COUNTRY_NAME}_{safe_id}_{year}.nc"
                if edd_outfile.exists():
                    yield edd_outfile
                continue

            print(f"  {item.id} ({year})")

            # --- Graceful skip on model failure ---
            try:
                tasmax = get_signed_dataset(item, "tasmax")["tasmax"]
                tasmin = get_signed_dataset(item, "tasmin")["tasmin"]
            except Exception as e:
                print(f"  ⚠️  {item.source_id} ({year}) — search/download failed: {e}, skipping")
                continue

            # lon fix + bbox slice
            tasmax = lon360_to_lon180(tasmax)
            tasmin = lon360_to_lon180(tasmin)
            ds_tasmax = tasmax.sel(region_bbox)
            ds_tasmin = tasmin.sel(region_bbox)

            print(f"    tasmax: {float(ds_tasmax.min()):.1f} – {float(ds_tasmax.max()):.1f} K")
            print(f"    tasmin: {float(ds_tasmin.min()):.1f} – {float(ds_tasmin.max()):.1f} K")

            if ds_tasmax.sizes.get("time", 0) == 0 or ds_tasmin.sizes.get("time", 0) == 0:
                print("    no data in bbox, skipping")
                continue

            # Kelvin → Celsius
            ds_tasmax = ds_tasmax - 273.15
            ds_tasmin = ds_tasmin - 273.15
            ds_tasmax.attrs["units"] = "Celsius"
            ds_tasmin.attrs["units"] = "Celsius"

            # Compute EDD
            try:
                edd_ds = calculate_daily_edd(ds_tasmax, ds_tasmin)
            except Exception as e:
                print(f"    EDD calculation failed: {e}, skipping")
                continue

            # Save EDD
            safe_id = str(item.id).replace("/", "_")
            edd_outfile = EDD_OUTPUT_DIR / f"edd_{COUNTRY_NAME}_{safe_id}_{year}.nc"

            try:
                edd_ds.to_netcdf(edd_outfile)
                print(f"    ✅ saved {edd_outfile.name}")
                save_checkpoint(item.source_id, year)   # mark as done
                yield edd_outfile
            except Exception as e:
                print(f"    ERROR saving: {e}")
                continue


# -------------------------------------------------------------------
# PART 2: Precompute calendar windows for this country
# -------------------------------------------------------------------

def compute_windows_for_calendar(gdf: gpd.GeoDataFrame, cal_name: str):
    dscal_path = CALENDAR_DIR / f"{cal_name}.crop.calendar.fill.nc"
    dscal = xr.open_dataset(str(dscal_path))

    gdf = gdf.reset_index(drop=True)   # ensure 0-based index
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

    T = 365
    offset = T

    for ii, row in gdf.iterrows():
        try:
            pix = dscal.sel(
                longitude=row["centroid"].x,
                latitude=row["centroid"].y,
                method="nearest",
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
    gdf = gpd.read_file(SHAPEFILE)
    gdf = gdf.to_crs("EPSG:4326")
    gdf = gdf.reset_index(drop=True)
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
