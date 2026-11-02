## config.R -- the only file you need to edit.
## DATA_DIR holds the SCImago files: sjr2004.csv ... sjr2025.csv and sjr_OA/sjr2004_OA.csv ... sjr2025_OA.csv
DATA_DIR <- "data"
OUT_DIR  <- file.path(DATA_DIR, "built")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
O <- function(f) file.path(OUT_DIR, f)
