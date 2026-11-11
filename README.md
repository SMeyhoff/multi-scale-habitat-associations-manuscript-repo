# Code and Supporting Information: **Spatiotemporal dynamics in scale of effect: multi-scale habitat associations of three large grazers across the growing season**

**Authors:** Sejer Meyhoff, Juan Andrés Martínez-Lanfranco, Isabelle Lebeuf-Taylor, Elly Knight, Erin Bayne (Department of Biological Sciences, University of Alberta; Elly Knight, Alberta Biodiversity Monitoring Institute). Contact: Sejer Meyhoff, meyhoff@ualberta.ca.

**Article:** accepted for publication in *Journal of Animal Ecology* (2026).

---

## Overview

This repository contains the full analytical workflow supporting the article. It is organized as eight self-contained R Markdown supporting-information documents (S1 to S8), three Jupyter notebooks that document the upstream PlanetScope NDVI processing pipeline, and the script that draws the methods schematic.

- **S1 to S6** are short, focused supplements that each address one methodological issue (temporal consistency of the NDVI cloud filter, NDVI and biomass validation, predictor screening and budgets, year-agnosticism of the model structure, the negative binomial choice, and spatial autocorrelation sensitivity).
- **S7** is the main analytical workflow: predictor construction, multivariate negative binomial GLMs across 16 buffer radii, scale-of-effect profiling, and the manuscript's result figures and tables.
- **S8** reports the sensitivity of scale-of-effect inference to specifying the vegetation stratum term as a fixed rather than a random effect. The fixed effect is the specification used throughout the manuscript.
- **Python/** contains two unrelated things. The three Jupyter notebooks are the NDVI pre-processing pipeline that produced the bi-weekly mosaics upstream of this analysis; they are included for transparency only, because the buffer-averaged NDVI values are already computed in `analysis_ready_multiscale.csv`. `methods_schematic.py` generates the article's methods schematic (Figure 2; see "Regenerating the methods schematic" below).

---

## Directory structure

```
repo/
├── README.md                              (this file)
├── LICENSE                                (MIT: applies to the code)
├── LICENSE.data                           (CC0: applies to the CSVs in input/)
├── S1_temporal_consistency_algorithm.Rmd  PlanetScope NDVI cloud filter: algorithm description
├── S2_ndvi_biomass_validation.Rmd         Exclosure NDVI and biomass validation
├── S3_predictor_budgets.Rmd               Predictor screening, budgets and selection
├── S4_year_agnostic_models.Rmd            Check that pooling years is appropriate
├── S5_overdispersion_nb_choice.Rmd        Justification for the negative binomial error family
├── S6_sac_sensitivity.Rmd                 Spatial autocorrelation sensitivity analysis
├── S7_Main_analysis.Rmd                   Main manuscript analysis and figures
├── S8_strata_specification.Rmd            Vegetation stratum as a fixed vs. random effect
│
├── Python/                                Non-R code
│   ├── phase1_8band_processing.ipynb      PlanetScope 8-band scene prep (transparency only)
│   ├── phase2_cloud_filtering.ipynb       Temporal-consistency cloud filter (see S1)
│   ├── monthly_ndvi_mosaics.ipynb         Bi-weekly and monthly mosaicing
│   └── methods_schematic.py               Generator for the methods schematic (Figure 2)
│
├── input/                                 Bundled data (CSVs; also deposited on Dryad)
│   ├── all_exclosure_data_by_year.csv     and 7 other CSVs, see "Data" below
│   ├── analysis_ready_multiscale.csv
│   ├── …                                  (8 CSVs total)
│   ├── cow.png                            Species icons used in the domain bar chart (Figure 3)
│   ├── elk.png
│   └── horse.png
│
├── figures/
│   └── range_exclosures.jpg               Photograph of a grazing exclosure, shown in S2
│
├── output/scale_analysis/                 CSV outputs written by S7 (domain classification,
│   │                                      top models by extent, stratum contrasts, and the
│   │                                      horse univariate comparison)
│   └── figures/                           Pre-rendered PNGs embedded in S7 (committed so the
│                                          HTML renders correctly on a fresh clone before a
│                                          full knit is run)
│
├── output/strata_specification/           CSV outputs written by S8 (coefficients, profiles,
│                                          variance components, VIFs and diagnostics under
│                                          both stratum specifications)
│
└── output/methods_schematic/              Written by Python/methods_schematic.py
    ├── methods_schematic.svg              Vector original
    └── methods_schematic.png              Raster version, as embedded in the manuscript
```

### Regenerating the methods schematic

The methods schematic is drawn programmatically rather than in a graphics editor, so it can be
regenerated from this repository. From `Python/`:

```
python methods_schematic.py ../output/methods_schematic/methods_schematic
```

Python 3.8+ standard library writes the SVG. `cairosvg` and `Pillow` are optional dependencies
used only to write the PNG.

The example scale-of-effect profile shown in the figure's result panel is **read from
`output/scale_analysis/domain_classification_full.csv`**, not transcribed, so the figure tracks the
fitted models. The script also holds an embedded copy of those values as a fallback and stops with
an error if the CSV disagrees with it, which prevents the figure drifting out of step with the
analysis after a re-knit. If the profile legitimately changes, update the embedded values in the
script after confirming the new knit is the one to publish.

---

## Reproducing the analysis

### 1. Data

The eight CSV files are bundled in `input/` so this repository is self-contained:

```
repo/input/
├── all_exclosure_data_by_year.csv        (read by S2)
├── analysis_ready_multiscale.csv         (read by S3, S4, S5, S6, S7, S8)
├── exclosure_ndvi_2023.csv               (read by S2)
├── exclosure_ndvi_2024.csv               (read by S2)
├── exclosure_sites.csv                   (read by S2)
├── exclosure_sites_ndvi.csv              (identical to exclosure_sites.csv; not read by any document, retained for legacy compatibility)
├── gt_deployment_effort.csv              (source records; not read by any document)
└── independent_detections_fhp_gt.csv     (source records; not read by any document)
```

`gt_deployment_effort.csv` and `independent_detections_fhp_gt.csv` are the upstream records from which the camera effort, the focal-species detection counts, and the human use predictor in `analysis_ready_multiscale.csv` were derived. The detections file is restricted to the three focal species and to the human and vehicle detections used to build the human use predictor. See the Dryad README for column definitions.

All Rmds use relative paths (`input/filename.csv`), so the files are picked up automatically with no working-directory configuration.

> **Canonical data citation.** The copies bundled here are provided for reproducibility convenience. The **authoritative, citable data source is the Dryad deposit** at DOI: 10.5061/dryad.vhhmgqp8k. If you reuse these data for other purposes, please cite the Dryad package, not this repository. The Dryad deposit also contains the full column-level data dictionary in its own README.
>
> **Licensing.** The CSVs in `input/` are released under **CC0 1.0 Universal** (see `LICENSE.data`), the same licence under which they are deposited on Dryad. The MIT licence in `LICENSE` applies only to the R and Python source code in this repository. Neither licence applies to the figures, the photograph, or the icon PNGs (see "License" below).

### 2. R environment

The analysis was run in **R 4.6.1** on Windows; the package versions used are listed in the session information at the end of S7 (`S7_Main_analysis.html`). The CRAN packages required to knit all eight Rmds:

```r
install.packages(c(
  "tidyverse",      # dplyr, tidyr, ggplot2, readr, purrr, stringr, lubridate, tibble
  "MASS",           # glm.nb
  "glmmTMB",        # negative binomial model fitting
  "DHARMa",         # residual diagnostics
  "MuMIn",          # dredge / model.sel
  "broom",
  "broom.mixed",
  "car",            # VIF / Anova
  "spdep",          # Moran's I, neighbour lists
  "furrr",          # parallel purrr
  "colorspace",
  "ggnewscale",
  "patchwork",
  "png",
  "ggstar",
  "kableExtra",
  "flextable",
  "officer",
  "knitr",
  "rmarkdown"
))
```

`parallel`, `grid` and `stats` ship with base R.

### 3. Knitting

From the `repo/` directory:

```r
rmarkdown::render("S7_Main_analysis.Rmd")
```

Caching is turned off in every document, so each knit refits all models from the data. S7 re-runs the full model set and took roughly ten minutes on the authors' workstation; S8 refits the base models under both stratum specifications. S1 to S6 are lighter.

Knit each document in its own R session (RStudio's Knit button does this). S4 loads `MASS`, whose `select()` masks `dplyr::select()` for any document rendered afterwards in the same session.

### 4. Pre-rendered HTML

If you would rather not re-knit, the rendered HTML for all eight documents is provided as Supporting Information with the manuscript.

---

## A note on the Python pipeline

The three notebooks in `Python/` document the upstream PlanetScope processing that produced the monthly NDVI rasters used to compute the buffer-averaged `NDVI_herb` values in `analysis_ready_multiscale.csv`. They require raw PlanetScope imagery (not provided: proprietary data licensed from Planet Labs) and a geospatial Python stack (Python 3.11) including `rasterio`, `numpy`, `pandas`, and `geopandas`.

The notebooks are included for **methodological transparency**: they let a reader inspect how the cloud filter described in S1 is implemented in code. They are **not part of the reproducibility workflow** for the article's models, because the buffer-averaged NDVI values are already computed and stored in the Dryad CSVs.

Absolute paths to the original network workspace have been replaced with a generic `<PROJECT_ROOT>/` token.

---

## License

- **Code** (the `.Rmd` files, the three Jupyter notebooks and the script in `Python/`) is released under the **MIT License**; see `LICENSE`.
- **Data** (the eight CSV files in `input/`) is released under **CC0 1.0 Universal (Public Domain Dedication)**; see `LICENSE.data`. This matches the licence under which the same data are deposited on Dryad at DOI: 10.5061/dryad.vhhmgqp8k.
- **Figures, photograph and icons** (`output/scale_analysis/figures/`, `output/methods_schematic/`, `figures/range_exclosures.jpg`, and `input/cow.png`, `input/elk.png`, `input/horse.png`) are provided so that the documents render and the analysis can be reproduced. They are not released under CC0 or MIT. Reuse terms for figures that appear in the article follow the licence of the published article.

---

## Citation

Please cite both the article and the data archive:

- **Article:** Meyhoff, S., Martínez-Lanfranco, J. A., Lebeuf-Taylor, I., Knight, E., & Bayne, E. (2026). Spatiotemporal dynamics in scale of effect: multi-scale habitat associations of three large grazers across the growing season. *Journal of Animal Ecology*.
- **Data:** Meyhoff, S., Martínez-Lanfranco, J. A., Lebeuf-Taylor, I., Knight, E., & Bayne, E. (2026). Data for: Spatiotemporal dynamics in scale of effect: multi-scale habitat associations of three large grazers across the growing season. Dryad Digital Repository. https://doi.org/10.5061/dryad.vhhmgqp8k
- **Code:** this repository, https://github.com/SMeyhoff/multi-scale-habitat-associations-manuscript-repo. Each release is archived on Zenodo; the Zenodo record gives the DOI to cite for the version used.
