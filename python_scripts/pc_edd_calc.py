import xarray as xr
import numpy as np

cutoff_temps = [0, 4, 8, 12, 28, 30, 32]

def calculate_daily_edd(ds_tasmax, ds_tasmin):
    """
    Calculates daily Exceedance Degree Days (EDD) using sinusoidal approximation.
    Based on Snyder (1985), Equation 4.

    Parameters:
        ds_tasmax: xarray.DataArray of daily max temperatures (°C)
        ds_tasmin: xarray.DataArray of daily min temperatures (°C)

    Returns:
        xarray.Dataset with EDDs for each threshold in cutoff_temps.
    """

    edd_results = {}

    # Midpoint and amplitude of the sine wave
    M = (ds_tasmax + ds_tasmin) / 2
    W = (ds_tasmax - ds_tasmin) / 2

    for t in cutoff_temps:
        R = (2 * t - ds_tasmax - ds_tasmin) / (ds_tasmax - ds_tasmin)
        R = R.clip(-1, 1)  # Ensure values are in valid domain for acos and sqrt

        theta = np.arccos(R)

        edd = xr.where(
            ds_tasmax <= t, 0,
            xr.where(
                ds_tasmin >= t, M - t,
                ((ds_tasmax - ds_tasmin) * theta / np.pi) -
                ((t - ds_tasmin) / np.pi) * np.sqrt(1 - R**2)
            )
        )

        edd_results[f"edd_{t}"] = edd

    edd_ds = xr.Dataset(edd_results)
    edd_ds.attrs["Conventions"] = "CF-1.6"
    edd_ds.attrs["title"] = "Daily Exceedance Degree-Days (Sinusoidal)"
    edd_ds.attrs["institution"] = "UD"
    edd_ds.attrs["source"] = "PC/CMIP6 tasmax + tasmin"

    return edd_ds


def calculate_monthly_edd(edd_ds):
    return edd_ds.resample(time="ME").mean()

