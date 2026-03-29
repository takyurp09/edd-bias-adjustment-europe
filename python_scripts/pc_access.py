import time
import planetary_computer
import pystac_client
import xarray as xr
import numpy as np


def get_items(source_ids, experiment_ids):
    print(f"  [PC] Searching catalog: source_ids={source_ids}, experiment_ids={experiment_ids}")
    catalog = pystac_client.Client.open(
        "https://planetarycomputer.microsoft.com/api/stac/v1",
        modifier=planetary_computer.sign_inplace,
    )
    search_results = catalog.search(
        collections=["cil-gdpcir-cc-by"],
        query={
            "cmip6:source_id": {"in": source_ids},
            "cmip6:experiment_id": {"in": experiment_ids},
        }
    )
    items = search_results.item_collection()
    print(f"  [PC] Found {len(items)} items")
    return items


def _to_lon180(ds: xr.Dataset) -> xr.Dataset:
    """Convert 0–360 longitude to -180–180."""
    if "lon" in ds.coords:
        lon = ds["lon"]
        if float(lon.max()) > 180:
            lon180 = ((lon + 180) % 360) - 180
            ds = ds.assign_coords(lon=lon180).sortby("lon")
    return ds


def _bbox_to_lon360(lon_min: float, lon_max: float):
    """
    Convert a -180/180 lon range to 0-360, handling the antimeridian.
    Returns a list of one or two (lon_min360, lon_max360) tuples.
    """
    lo = lon_min % 360
    hi = lon_max % 360
    if lo <= hi:
        return [(lo, hi)]
    else:
        return [(lo, 360.0), (0.0, hi)]


def get_signed_dataset(item, asset_key, year=None, region_bbox=None):
    """
    Opens a Planetary Computer Zarr asset efficiently by:
      1. Subsetting lon/lat in native 0-360 space (chunk-aligned)
      2. Subsetting time before any compute
      3. Converting lon to -180/180 only on the already-small array
    """
    print(f"    [PC] Signing asset: {asset_key}")
    t_sign = time.time()
    signed_asset = planetary_computer.sign(item.assets[asset_key])
    print(f"    [PC] Signed in {time.time()-t_sign:.1f}s")

    open_kwargs = signed_asset.extra_fields.get("xarray:open_kwargs", {})
    href = signed_asset.href

    print(f"    [PC] Opening Zarr store...")
    t_open = time.time()
    try:
        ds = xr.open_zarr(href, consolidated=True, **open_kwargs)
    except Exception:
        print(f"    [PC] open_zarr failed, falling back to open_dataset")
        ds = xr.open_dataset(href, **open_kwargs)
    print(f"    [PC] Zarr store opened in {time.time()-t_open:.1f}s")

    full_time = ds.sizes.get("time", "?")
    full_lon  = ds.sizes.get("lon",  "?")
    full_lat  = ds.sizes.get("lat",  "?")
    print(f"    [PC] Full store shape — time:{full_time}  lat:{full_lat}  lon:{full_lon}")

    # --- Time subset ---
    if year is not None:
        ds = ds.sel(time=slice(f"{year}-01-01", f"{year}-12-31"))
        print(f"    [PC] After time slice ({year}): {ds.sizes.get('time', '?')} timesteps")

    # --- Spatial subset in NATIVE lon coords (0-360) ---
    if region_bbox is not None:
        lat_slice = region_bbox.get("lat")
        lon_slice = region_bbox.get("lon")

        if lat_slice is not None:
            ds = ds.sel(lat=lat_slice)

        if lon_slice is not None and "lon" in ds.coords:
            native_lon = ds["lon"]
            if float(native_lon.max()) > 180:
                lo180 = lon_slice.start
                hi180 = lon_slice.stop
                slices360 = _bbox_to_lon360(lo180, hi180)
                print(f"    [PC] native lon is 0-360, converting bbox to 360 space: {slices360}")
                if len(slices360) == 1:
                    lo360, hi360 = slices360[0]
                    ds = ds.sel(lon=slice(lo360, hi360))
                else:
                    parts = [ds.sel(lon=slice(lo, hi)) for lo, hi in slices360]
                    ds = xr.concat(parts, dim="lon")
            else:
                print(f"    [PC] native lon is already -180/180, subsetting directly")
                ds = ds.sel(lon=lon_slice)

        print(f"    [PC] After bbox subset: lat:{ds.sizes.get('lat','?')}  lon:{ds.sizes.get('lon','?')}")

    # --- Convert lon to -180/180 on the small regional array ---
    t_lon = time.time()
    ds = _to_lon180(ds)
    print(f"    [PC] lon conversion done in {time.time()-t_lon:.1f}s")
    print(f"    [PC] Ready — final shape: {dict(ds.sizes)}")

    return ds