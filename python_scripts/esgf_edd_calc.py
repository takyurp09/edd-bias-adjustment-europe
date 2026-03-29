import xarray as xr
import numpy as np
import shutil
from pathlib import Path
from pangeo_access import load_cmip6_year

# ----------------------------
# CONFIG
# ----------------------------
cutoff_temps = [0, 4, 8, 12, 28, 30, 32]
TMP_DIR = "data/raw_cmip6/tmp"

# ----------------------------
# CLEANUP
# ----------------------------
def cleanup_raw_cmip6(tmp_dir=TMP_DIR):
    p = Path(tmp_dir)
    if p.exists():
        shutil.rmtree(p)

# ----------------------------
# DAILY EDD (SINUSOIDAL)
# ----------------------------
def calculate_daily_edd(ds_tasmax, ds_tasmin):
    """
    Sinusoidal exceedance degree-days.
    IMPORTANT: inputs MUST already be in Celsius.
    """

    edd_results = {}

    # Mean and half-range
    M = (ds_tasmax + ds_tasmin) / 2.0
    W = (ds_tasmax - ds_tasmin) / 2.0

    for t in cutoff_temps:
        # Ratio for arccos, clipped for numerical safety
        R = (t - M) / W
        R = R.clip(-1.0, 1.0)

        theta = np.arccos(R)

        # Correct sinusoidal exceedance integral
        edd_mid = (W * np.sin(theta) + (M - t) * theta) / np.pi

        edd = xr.where(
            t >= ds_tasmax, 0.0,
            xr.where(
                t <= ds_tasmin, M - t,
                edd_mid
            )
        )

        # Handle flat days safely
        edd = xr.where(W == 0, xr.where(M > t, M - t, 0.0), edd)

        # Numerical safety
        edd = edd.clip(min=0.0)

        edd_results[f"edd_{t}"] = edd.astype("float32")

    edd_ds = xr.Dataset(edd_results)
    edd_ds.attrs["title"] = "Daily Exceedance Degree-Days (Sinusoidal, Celsius)"
    edd_ds.attrs["source"] = "CMIP6 tasmax/tasmin (ESGF)"
    edd_ds.attrs["institution"] = "University of Delaware"

    return edd_ds

# ----------------------------
# MONTHLY AGGREGATION
# ----------------------------
def calculate_monthly_edd(edd_ds):
    return edd_ds.resample(time="ME").sum()

# ----------------------------
# YEARLY WRAPPER (optional utility)
# ----------------------------
def run_edd_for_year(year, region_bbox=None):
    ds_tmax = load_cmip6_year(
        source_id="GFDL-ESM4",
        experiment_id="historical",
        table_id="day",
        variable_id="tasmax",
        year=year,
        region_bbox=region_bbox,
    )

    ds_tmin = load_cmip6_year(
        source_id="GFDL-ESM4",
        experiment_id="historical",
        table_id="day",
        variable_id="tasmin",
        year=year,
        region_bbox=region_bbox,
    )

    # Convert ONCE here
    tmax_c = ds_tmax["tasmax"] - 273.15
    tmin_c = ds_tmin["tasmin"] - 273.15

    edd_daily = calculate_daily_edd(tmax_c, tmin_c)
    edd_monthly = calculate_monthly_edd(edd_daily)

    cleanup_raw_cmip6()
    return edd_daily, edd_monthly
