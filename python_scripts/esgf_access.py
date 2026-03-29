import os
import time
import intake_esgf
import xarray as xr
from intake_esgf.exceptions import NoSearchResults

# ---------- helpers (unchanged) ----------

class ESGFItem:
    def __init__(self, source_id, experiment_id):
        self.source_id = source_id
        self.experiment_id = experiment_id
        self.id = f"{source_id}_{experiment_id}"
        self.assets = {"tasmax": True, "tasmin": True, "pr": True}

def get_items(source_ids, experiment_ids):
    items = []
    for s in source_ids:
        for e in experiment_ids:
            items.append(ESGFItem(s, e))
    return items

def _to_lon180(ds: xr.Dataset) -> xr.Dataset:
    if "lon" in ds.coords:
        lon = ds["lon"]
        if float(lon.max()) > 180:
            lon180 = ((lon + 180) % 360) - 180
            ds = ds.assign_coords(lon=lon180).sortby("lon")
    return ds

def _pick_grid(df):
    """Prefer native grid (gn), then regridded (gr/gr1), then first available."""
    if "grid_label" in df.columns:
        for label in ("gn", "gr1", "gr"):
            sub = df[df["grid_label"] == label]
            if not sub.empty:
                return sub.iloc[[0]]
    return df.iloc[[0]]

# ---------- search strategies ----------

def _search_strict(cat, source_id, experiment_id, variable_id, member_id):
    """Tight search: explicit member_id, no ignore_facets."""
    return cat.search(
        project="CMIP6",
        source_id=source_id,
        experiment_id=experiment_id,
        table_id="day",
        variable_id=variable_id,
        member_id=member_id,
    )

def _search_any_member(cat, source_id, experiment_id, variable_id):
    """Relax member_id constraint — accept any ensemble member."""
    return cat.search(
        project="CMIP6",
        source_id=source_id,
        experiment_id=experiment_id,
        table_id="day",
        variable_id=variable_id,
    )

def _search_ignore_facets(cat, source_id, experiment_id, variable_id):
    """Last resort: ignore_facets=True (original behaviour)."""
    return cat.search(
        project="CMIP6",
        source_id=source_id,
        experiment_id=experiment_id,
        table_id="day",
        variable_id=variable_id,
        ignore_facets=True,
    )

# ---------- main public functions ----------

def get_signed_dataset(
    item,
    asset_key,
    region_bbox=None,
    member_id: str = "r1i1p1f1",
    max_retries: int = 3,
    retry_delay: float = 10.0,
):
    year = int(os.environ["CMIP6_YEAR"])

    strategies = [
        ("strict",        lambda c: _search_strict(c, item.source_id, item.experiment_id, asset_key, member_id)),
        ("any_member",    lambda c: _search_any_member(c, item.source_id, item.experiment_id, asset_key)),
        ("ignore_facets", lambda c: _search_ignore_facets(c, item.source_id, item.experiment_id, asset_key)),
    ]

    last_exc = None

    for attempt in range(1, max_retries + 1):
        for strategy_name, search_fn in strategies:
            try:
                cat = intake_esgf.ESGFCatalog()
                search = search_fn(cat)

                df = search.df
                if df is None or df.empty:
                    print(f"  [{strategy_name}] empty df, skipping")
                    continue

                search.df = _pick_grid(df)
                dsets = search.to_dataset_dict()

                if not dsets:
                    print(f"  [{strategy_name}] to_dataset_dict returned nothing, skipping")
                    continue

                key = list(dsets.keys())[0]
                ds = dsets[key]
                print(f"  [{strategy_name}] SUCCESS — source key: {key}")

                ds = ds.sel(time=slice(f"{year}-01-01", f"{year}-12-31"))
                ds = _to_lon180(ds)
                if region_bbox is not None:
                    ds = ds.sel(region_bbox)

                return ds

            except NoSearchResults:
                print(f"  [{strategy_name}] NoSearchResults, trying next strategy")
                last_exc = NoSearchResults()
            except Exception as exc:
                print(f"  [{strategy_name}] unexpected error: {exc}")
                last_exc = exc

        # All strategies failed this attempt — wait before retrying
        if attempt < max_retries:
            print(f"  All strategies failed (attempt {attempt}/{max_retries}). "
                  f"Retrying in {retry_delay}s …")
            time.sleep(retry_delay)

    raise ValueError(
        f"No dataset found for {item.id} / {asset_key} after {max_retries} attempts. "
        f"Last error: {last_exc}"
    )


def load_cmip6_year(
    source_id: str,
    experiment_id: str,
    variable_id: str,
    year: int,
    table_id: str = "day",
    member_id: str = "r1i1p1f1",
    project: str = "CMIP6",
    region_bbox=None,
    max_retries: int = 3,
    retry_delay: float = 10.0,
) -> xr.Dataset:

    strategies = [
        ("strict",     lambda c: c.search(project=project, source_id=source_id,
                                           experiment_id=experiment_id, table_id=table_id,
                                           variable_id=variable_id, member_id=member_id)),
        ("any_member", lambda c: c.search(project=project, source_id=source_id,
                                           experiment_id=experiment_id, table_id=table_id,
                                           variable_id=variable_id)),
    ]

    last_exc = None

    for attempt in range(1, max_retries + 1):
        for strategy_name, search_fn in strategies:
            try:
                cat = intake_esgf.ESGFCatalog()
                search = search_fn(cat)

                df = search.df
                if df is None or df.empty:
                    continue

                search.df = _pick_grid(df)
                dsets = search.to_dataset_dict(aggregate=False)

                if not dsets:
                    continue

                ds = next(iter(dsets.values()))
                print(f"  [{strategy_name}] SUCCESS — {list(dsets.keys())[0]}")

                ds = ds.sel(time=slice(f"{year}-01-01", f"{year}-12-31"))
                ds = _to_lon180(ds)
                if region_bbox is not None:
                    ds = ds.sel(region_bbox)

                return ds

            except NoSearchResults:
                print(f"  [{strategy_name}] NoSearchResults, trying next strategy")
                last_exc = NoSearchResults()
            except Exception as exc:
                print(f"  [{strategy_name}] unexpected error: {exc}")
                last_exc = exc

        if attempt < max_retries:
            print(f"  All strategies failed (attempt {attempt}/{max_retries}). "
                  f"Retrying in {retry_delay}s …")
            time.sleep(retry_delay)

    raise ValueError(
        f"No dataset found for {source_id} {experiment_id} {variable_id} "
        f"after {max_retries} attempts. Last error: {last_exc}"
    )