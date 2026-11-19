# How Closely Does SCImago Journal Rank Track Citations?

Code for: Kim, E. *How Closely Does SCImago Journal Rank Track Citations? Within-Journal Citation Elasticity Across Fields, Publishers, and Years*.

## Data

Download the yearly files from [SCImago Journal & Country Rank](https://www.scimagojr.com/journalrank.php) and put them in `data/`:

```
data/sjr2004.csv ... data/sjr2025.csv
```

## Run

```
Rscript 1_build_panel.R      # panel and samples
Rscript 2_models.R           # Models M1–M6, scalings, robustness
Rscript 3_scaling_checks.R   # conversion to elasticities, publisher gap
python  4_fe_checks.py       # field × year fixed effects, composition and journal-size checks
Rscript 5_figures.R          # Figures 1–4
```

Results go to `data/built/`, figures to `figures/`. The full run takes about 20–30 minutes.

## Software

R 4.3 (data.table, stringr, plm, lmtest, sandwich, jsonlite, ggplot2, patchwork, ggrepel); Python 3 (pandas, numpy, pyfixest 0.60).

## License

MIT. SCImago data are subject to SCImago's terms of use.
