# What Does SCImago Journal Rank Reward?

Code for: Kim, E. *What Does SCImago Journal Rank Reward? Within-Journal Citation Elasticity Across Fields, Publishers, and Years, 2004–2025*. Research Evaluation (under review).

## Data

Download the yearly files from [SCImago Journal & Country Rank](https://www.scimagojr.com/journalrank.php) and put them in `data/`:

```
data/sjr2004.csv ... data/sjr2025.csv                 all journals
data/sjr_OA/sjr2004_OA.csv ... sjr_OA/sjr2025_OA.csv  "Only Open Access Journals" filter (optional)
```

Or copy the `built/` folder of the data deposit to `data/built/` and start at step 2.

## Run

```
Rscript 1_build_panel.R      # panel and samples
Rscript 2_models.R           # Models M1–M8, scalings, robustness
Rscript 3_scaling_checks.R   # conversion to elasticities, publisher gap
python  4_fe_checks.py       # field × year fixed effects, composition checks
Rscript 5_figures.R          # Figures 1–4
```

Results go to `data/built/`, figures to `figures/`. The full run takes about 20–30 minutes.

## Software

R 4.3 (data.table, stringr, plm, lmtest, sandwich, jsonlite, ggplot2, patchwork, ggrepel); Python 3 (pandas, numpy, pyfixest 0.60).

## License

MIT. SCImago data are subject to SCImago's terms of use.
