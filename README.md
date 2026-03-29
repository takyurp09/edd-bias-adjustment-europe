# Evaluating Bias-Adjusted Climate Data for Extreme Heat Metrics in Agricultural Applications

**Paper**: Tahmid, M. T. (2026). Evaluating Bias-Adjusted Climate Data for Extreme Heat Metrics in Agricultural Applications. *Environmental Data Science*.

## Overview

This repository contains code for reproducing the analysis comparing raw CMIP6, bias-adjusted CIL-GDPCIR, and ERA5 reanalysis for Extreme Degree Days (EDD) metrics across European agriculture (1994-2014).

## Repository Structure
```
.
├── README.md                          # This file
├── LICENSE                            # MIT License
├── environment.yml                    # Conda environment
├── 00_config.R                        # Configuration
├── 01_load_harmonize.R                # Load data
├── 02_metrics.R                       # Calculate metrics
├── 03_maps.R                          # Generate maps
├── 04_plots.R                         # Generate plots
├── 05_generate_paper_tables.R         # Extract statistics
├── 06_generate_conference_figures.R   # Generate figures
├── run_all.R                          # Master script
└── python_scripts/                    # Python code
    ├── era5_*.py                      # ERA5 processing
    ├── esgf_*.py                      # Raw CMIP6 access
    └── pc_*.py                        # Bias-adjusted CMIP6 access
```

## Requirements

- Python 3.9+
- R 4.2+
- See `environment.yml` for dependencies

## Installation
```bash
conda env create -f environment.yml
conda activate edd-analysis
```

## Data

**Aggregated results** (~500 MB) available separately on Zenodo: [DOI will be added]

**Raw climate data** (not included - cite sources):
- ERA5: Copernicus Climate Data Store
- Raw CMIP6: Pangeo archive
- Bias-adjusted CMIP6: Microsoft Planetary Computer (CIL-GDPCIR)
- Crop calendars: SAGE (Sacks et al. 2010)
- Harvested area: GAEZ v4.0 (Fischer et al. 2021)

## Citation
```bibtex
@article{tahmid2026,
  author = {Tahmid, Muhammad Taky},
  title = {Evaluating Bias-Adjusted Climate Data for Extreme Heat Metrics in Agricultural Applications},
  journal = {Environmental Data Science},
  year = {2026}
}
```

## License

MIT License

## Contact

Muhammad Taky Tahmid - tahmid@udel.edu
