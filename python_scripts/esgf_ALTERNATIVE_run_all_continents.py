#!/usr/bin/env python

import os
import subprocess

# (name, tag, adm_level)
"""
CONTINENTS = [
    ("africa",        "AFR", 1),
    ("asia",          "ASIA", 1),
    ("europe",        "EUR", 1),
    ("north_america", "NAM", 1),
    ("south_america", "SAM", 1),
    ("oceania",       "OCE", 1),
]
"""

CONTINENTS = [
    ("europe",        "EUR", 1)
]



for country_name, country_tag, adm_level in CONTINENTS:
    print("\n" + "=" * 70)
    print(f"Running for {country_name} ({country_tag}), ADM_LEVEL={adm_level}")
    print("=" * 70 + "\n")

    env = os.environ.copy()
    env["COUNTRY_NAME"] = country_name
    env["COUNTRY_TAG"] = country_tag
    env["ADM_LEVEL"] = str(adm_level)

    # If you normally run `python ALTERNATIVE_agg_edd.py`, keep that here
    result = subprocess.run(
        ["python", "esgf_ALTERNATIVE_agg_edd.py"],
        env=env,
        check=False,
    )

    if result.returncode != 0:
        print(f"\n❌ Stopping: {country_name} failed with return code {result.returncode}")
        break

    print(f"\n✅ Finished {country_name} ({country_tag})")
