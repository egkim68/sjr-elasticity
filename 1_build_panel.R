## =============================================================
## 1_build_panel.R
## Builds the journal-year panel from the 22 annual SCImago files,
## applies the sample filters, and writes the analysis files to data/built/.
## =============================================================
suppressMessages({library(data.table); library(stringr)})
source("config.R")
YEARS <- 2004:2025

## ---------- helpers
split_issns <- function(x) {
  x <- str_split(x, "[,;]")[[1]]
  x <- str_trim(x)
  x <- str_replace_all(x, "-", "")
  x[x != ""]
}

norm_title <- function(x) {
  x <- tolower(x)
  x <- str_replace_all(x, "[^a-z0-9]+", " ")
  str_trim(x)
}

read_sjr_year <- function(path, year) {
  dt <- fread(
    path,
    sep = ";",
    encoding = "UTF-8",
    colClasses = "character",
    na.strings = c("", "NA")
  )
  ## the export has two "Publisher" columns; the second becomes "Publisher.1"
  setnames(dt, make.unique(names(dt)))
  dt[, year := year]
  dt
}

## "145,004" -> 145.004
sjr_num <- function(x) as.numeric(str_replace_all(x, ",", "."))

## ---------- read the 22 annual files
main_list <- vector("list", length(YEARS))
names(main_list) <- YEARS

for (y in YEARS) {
  f <- file.path(DATA_DIR, sprintf("sjr%d.csv", y))
  if (!file.exists(f)) {
    warning(sprintf("Main file missing for %d, skipping: %s", y, f))
    next
  }
  dt <- read_sjr_year(f, y)
  ## "Total Docs. (2024)" -> "Total Docs. (year)"
  docs_col <- grep("^Total Docs\\. \\(\\d{4}\\)$", names(dt), value = TRUE)
  if (length(docs_col) == 1) setnames(dt, docs_col, "Total Docs. (year)")

  main_list[[as.character(y)]] <- dt
}

main_all <- rbindlist(main_list, use.names = TRUE, fill = TRUE)
rm(main_list)

## ---------- read the open-access lists
oa_list <- vector("list", length(YEARS))
names(oa_list) <- YEARS

for (y in YEARS) {
  f <- file.path(DATA_DIR, "sjr_OA", sprintf("sjr%d_OA.csv", y))
  if (!file.exists(f)) {
    warning(sprintf("OA file missing for %d, skipping: %s", y, f))
    next
  }
  dt <- read_sjr_year(f, y)
  oa_list[[as.character(y)]] <- dt
}

oa_all <- rbindlist(oa_list, use.names = TRUE, fill = TRUE)
rm(oa_list)
## without the sjr_OA/ folder, the native "Open Access" flag is used instead
OA_FILES_PRESENT <- nrow(oa_all) > 0

## ---------- OA status: match each year's OA list by ISSN, then by title
build_issn_lookup <- function(dt) {
  out <- dt[, .(issn = split_issns(Issn)), by = .(year, Sourceid, Title)]
  out
}

if (OA_FILES_PRESENT) {
main_issn_long <- build_issn_lookup(main_all)
oa_issn_long   <- build_issn_lookup(oa_all)

oa_issn_key <- unique(oa_issn_long[, .(year, issn)])
oa_issn_key[, oa_by_issn := TRUE]

main_issn_long <- merge(main_issn_long, oa_issn_key,
                         by = c("year", "issn"), all.x = TRUE)

matched_by_issn <- unique(
  main_issn_long[oa_by_issn == TRUE, .(year, Sourceid)]
)
matched_by_issn[, oa_by_issn := TRUE]

main_all[, title_norm := norm_title(Title)]
oa_all[,   title_norm := norm_title(Title)]

oa_title_key <- unique(oa_all[, .(year, title_norm)])
oa_title_key[, oa_by_title := TRUE]

main_all <- merge(main_all, oa_title_key,
                   by = c("year", "title_norm"), all.x = TRUE)
main_all <- merge(main_all, matched_by_issn,
                   by = c("year", "Sourceid"), all.x = TRUE)

main_all[, oa_status := fifelse(
  !is.na(oa_by_issn) | !is.na(oa_by_title), 1L, 0L
)]
main_all[, c("oa_by_issn", "oa_by_title", "title_norm") := NULL]
} else {
  cat("sjr_OA/ files not found: using the native 'Open Access' flag as oa_status\n")
  main_all[, oa_status := fifelse(`Open Access` == "Yes", 1L, 0L)]
}

## agreement with the native "Open Access" flag, by year file (reported in Section 2.2)
main_all[, oa_native := fifelse(`Open Access` == "Yes", 1L, 0L)]
print(main_all[, .(agree_rate = mean(oa_status == oa_native, na.rm = TRUE)), by = year][order(year)])

## ---------- short variable names and numeric conversion
old_names <- c(
  "Rank", "Sourceid", "Title", "Type", "Issn", "Publisher",
  "Open Access Diamond", "SJR", "SJR Best Quartile", "H index",
  "Total Docs. (year)", "Total Docs. (3years)", "Total Refs.",
  "Total Citations (3years)", "Citable Docs. (3years)",
  "Citations / Doc. (2years)", "Ref. / Doc.", "%Female", "Overton",
  "Country", "Region", "Coverage", "Categories", "Areas"
)
new_names <- c(
  "rank", "src_id", "title", "src_type", "issn", "publ",
  "oa_diamond", "sjr", "quart", "hindex",
  "docs", "docs3y", "refs",
  "cites3y", "citdocs3y",
  "cpd", "rpd", "pctfem", "overton",
  "country", "region", "coverage", "categories", "areas"
)
keep <- old_names %in% names(main_all)
setnames(main_all, old_names[keep], new_names[keep])

num_cols <- c("sjr", "hindex", "docs", "docs3y", "refs", "cites3y",
              "citdocs3y", "cpd", "rpd", "pctfem", "overton")
num_cols <- intersect(num_cols, names(main_all))
main_all[, (num_cols) := lapply(.SD, sjr_num), .SDcols = num_cols]

## ---------- NYears: indexed years up to each panel year, from the Coverage ranges
## (e.g. "1973-1975, 1978-2026"; gaps are not counted, and ranges are clipped at the panel year)
cov_lookup <- unique(
  main_all[!is.na(coverage) & coverage != "", .(src_id, coverage)],
  by = "src_id"
)

range_pieces <- cov_lookup[, .(piece = str_trim(str_split(coverage, ",")[[1]])),
                            by = src_id]
range_pieces[, has_dash := str_detect(piece, "-")]
range_pieces[, `:=`(
  start = fifelse(has_dash,
                   as.integer(str_extract(piece, "^\\d{4}")),
                   as.integer(piece)),
  end   = fifelse(has_dash,
                   as.integer(str_extract(piece, "\\d{4}$")),
                   as.integer(piece))
)]
range_pieces <- range_pieces[!is.na(start) & !is.na(end), .(src_id, start, end)]
range_pieces[, start_orig := start]   # kept, since the join below overwrites start

panel_key <- unique(main_all[, .(src_id, year)])

joined <- range_pieces[panel_key, on = .(src_id, start <= year),
                        allow.cartesian = TRUE, nomatch = NULL]
joined[, contrib := pmin(end, start) - start_orig + 1L]
nyears_lookup <- joined[, .(nyears_covered = sum(contrib)), by = .(src_id, year = start)]

main_all <- merge(main_all, nyears_lookup, by = c("src_id", "year"), all.x = TRUE)

## ---------- born year from Coverage; journals covered from 1996 or earlier are "early journals"
born_year <- range_pieces[, .(born_year = min(start_orig)), by = src_id]
main_all <- merge(main_all, born_year, by = "src_id", all.x = TRUE)
main_all[, retro_flag := born_year <= 1996]
saveRDS(main_all, O("sjr_panel_full.rds"))

## ---------- sample: journals and trade journals with a valid quartile, early journals excluded
panel <- main_all[src_type %in% c("journal", "trade journal") & !is.na(quart) & !quart %in% c("", "-")]
cat("Journals with a valid quartile:", uniqueN(panel$src_id), "| early journals:",
    panel[retro_flag %in% TRUE, uniqueN(src_id)], "| blank coverage:", panel[is.na(retro_flag), uniqueN(src_id)], "\n")
panel <- panel[retro_flag == FALSE]
setorder(panel, src_id, year)
cat("Analytic sample:", uniqueN(panel$src_id), "journals,", nrow(panel), "journal-years\n")

## ---------- quartile persistence (share keeping the same quartile from one year to the next)
panel[, `:=`(quart_lag = shift(quart), year_lag = shift(year)), by = src_id]
trans <- panel[!is.na(quart_lag) & year == year_lag + 1]
fwrite(trans[, .(persistence_rate = mean(quart == quart_lag)), by = quart_lag][order(quart_lag)],
       O("table4a_persistence_by_quartile.csv"))
panel[, c("quart_lag", "year_lag") := NULL]

## ---------- analytic panel: log outcome and predictors lagged one year
panel[, log_sjr := log(sjr + 0.001)]
setorder(panel, src_id, year)
panel[, `:=`(
  docs_l1       = shift(docs, 1),
  nyears_l1     = shift(nyears_covered, 1),
  cpd_l1        = shift(cpd, 1),
  hindex_l1     = shift(hindex, 1),
  rpd_l1        = shift(rpd, 1),
  log_sjr_l1    = shift(log_sjr, 1),
  year_l1       = shift(year, 1)
), by = src_id]

## ---------- 27 field dummies (a journal can belong to several fields; exact match on area names)
all_fields <- sort(unique(str_trim(unlist(str_split(panel$areas, ";")))))
field_var_names <- paste0("field_", make.names(all_fields))
area_tokens <- lapply(str_split(panel$areas, ";"), str_trim)
for (i in seq_along(all_fields)) {
  panel[, (field_var_names[i]) := as.integer(vapply(area_tokens, function(t) all_fields[i] %in% t, logical(1)))]
}
rm(area_tokens)

## ---------- publishers: consolidate name variants and imprints, then keep those
## with at least 100 journals (21 publishers; publisher is a current snapshot)
panel[, Publisher.1 := str_replace_all(Publisher.1, fixed("&amp;"), "&")]
panel[, pub_group := Publisher.1]
panel[Publisher.1 %in% c("Elsevier B.V.","Elsevier Ltd","Elsevier Inc.","Elsevier GmbH",
  "Elsevier Masson s.r.l.","Elsevier Ireland Ltd","Elsevier (Singapore) Pte Ltd",
  "Elsevier Editora Ltda","Elsevier Espana S.L.U","Elsevier","Reed Elsevier India Pvt. Ltd.",
  "Academic Press Inc.","Academic Press","W.B. Saunders","W.B. Saunders Ltd",
  "Churchill Livingstone","Cell Press"), pub_group := "Elsevier"]
panel[Publisher.1 %in% c("Springer Science and Business Media Deutschland GmbH","Springer",
  "Springer Verlag","Springer International Publishing AG","Springer Netherlands",
  "Springer Science and Business Media B.V.","Springer New York","Springer London",
  "Springer Science + Business Media","Springer International Publishing","SpringerOpen",
  "Springer India","Springer Singapore","Springer Paris","Springer Heidelberg","Springer Japan",
  "Springer-Verlag Wien","Springer-Verlag Italia s.r.l.","Springer Vieweg","Springer US",
  "Springer Medizin","Springer GmbH & Co, Auslieferungs-Gesellschaf",
  "Springer Science and Business Media, LLC","Springer Nature"), pub_group := "Springer"]
panel[Publisher.1 %in% c("Nature Research","Nature Publishing Group"), pub_group := "Nature Research"]
panel[Publisher.1 == "BioMed Central Ltd", pub_group := "BioMed Central"]
panel[Publisher.1 == "Palgrave Macmillan Ltd.", pub_group := "Palgrave Macmillan"]
panel[Publisher.1 %in% c("John Wiley and Sons Inc","John Wiley and Sons Ltd",
  "John Wiley & Sons Inc.","John Wiley &amp; Sons Inc.","Wiley-Liss Inc."), pub_group := "Wiley"]
panel[Publisher.1 %in% c("Wiley-Blackwell Publishing Ltd","Wiley-Blackwell",
  "Blackwell Publishing Asia","Blackwell Munksgaard"), pub_group := "Wiley-Blackwell"]
panel[Publisher.1 %in% c("Wiley-VCH Verlag","Wiley - VCH Verlag GmbH & CO. KGaA",
  "Wiley - VCH Verlag GmbH &amp; CO. KGaA"), pub_group := "Wiley-VCH"]
panel[Publisher.1 %in% c("Taylor and Francis Ltd.","Taylor and Francis Inc.",
  "Taylor & Francis Group LLC","Taylor & Francis Group LLC Philadelphia",
  "Informa UK Ltd"), pub_group := "Taylor & Francis"]
panel[Publisher.1 %in% c("SAGE Publications Inc.","SAGE Publications Ltd",
  "Sage Publications India Pvt. Ltd","Sage Publications","Sage Periodicals Press"), pub_group := "SAGE Publications"]
panel[Publisher.1 %in% c("Lippincott Williams and Wilkins","Lippincott Williams and Wilkins Ltd."), pub_group := "Lippincott Williams & Wilkins"]
panel[Publisher.1 %in% c("Wolters Kluwer Medknow Publications","Wolters Kluwer Health",
  "Wolters Kluwer Health Inc"), pub_group := "Wolters Kluwer"]
panel[Publisher.1 %in% c("Institute of Electrical and Electronics Engineers Inc.",
  "Institute of Electrical and Electronics Engineers","IEEE Computer Society",
  "IEEE Computational Intelligence Society","IEEE Advancing Technology for Humanity",
  "IEEE Communications Society"), pub_group := "IEEE"]
panel[Publisher.1 %in% c("Walter de Gruyter GmbH","Walter de Gruyter GmbH & Co. KG",
  "De Gruyter Mouton","De Gruyter Open Ltd","Sciendo"), pub_group := "De Gruyter"]
panel[Publisher.1 %in% c("Emerald Group Publishing Ltd.","Emerald Publishing",
  "Emarald Group Publishing Ltd"), pub_group := "Emerald"]
panel[Publisher.1 %in% c("Multidisciplinary Digital Publishing Institute (MDPI)","MDPI AG"), pub_group := "MDPI"]
panel[Publisher.1 %in% c("Inderscience Enterprises Ltd","Inderscience","Inderscience Publishers"), pub_group := "Inderscience"]
panel[Publisher.1 %in% c("World Scientific","World Scientific Publishing Co. Pte Ltd"), pub_group := "World Scientific"]
panel[Publisher.1 %in% c("Akademiai Kiado","Akademiai Kiado ZRt."), pub_group := "Akademiai Kiado"]
panel[Publisher.1 %in% c("Kluwer Academic Publishers","Kluwer Academic/Human Sciences Press Inc."), pub_group := "Kluwer Academic Publishers"]
panel[Publisher.1 == "Kluwer Law International", pub_group := "Kluwer Law International"]
panel[Publisher.1 %in% c("Optica Publishing Group (formerly OSA)","Optica Publishing Group"), pub_group := "Optica Publishing Group"]
panel[Publisher.1 %in% c("KeAi Communications Co.","KeAi Publishing Communications Ltd."), pub_group := "KeAi"]
panel[Publisher.1 %in% c("Georg Thieme Verlag","Thieme Medical Publishers, Inc."), pub_group := "Thieme"]
panel[Publisher.1 %in% c("Fabrizio Serra Editore Srl","Fabrizio Serra Editore"), pub_group := "Fabrizio Serra Editore"]
panel[str_detect(Publisher.1, regex("Bentham Science", ignore_case = TRUE)), pub_group := "Bentham Science"]
panel[str_detect(Publisher.1, regex("Brill", ignore_case = TRUE)), pub_group := "Brill"]
panel[str_detect(Publisher.1, regex("Frontiers Media", ignore_case = TRUE)), pub_group := "Frontiers"]
panel[str_detect(Publisher.1, regex("Oxford University Press", ignore_case = TRUE)), pub_group := "Oxford University Press"]
panel[str_detect(Publisher.1, regex("Cambridge University Press", ignore_case = TRUE)), pub_group := "Cambridge University Press"]

pub_vol <- panel[pub_group != "", .(n_journals = uniqueN(src_id)), by = pub_group]
setorder(pub_vol, -n_journals)
publisher_groups <- pub_vol[n_journals >= 100]$pub_group
pub_var_names <- paste0("pub_", make.names(publisher_groups))
for (i in seq_along(publisher_groups)) {
  panel[, (pub_var_names[i]) := as.integer(pub_group == publisher_groups[i])]
}

## ---------- estimation sample: consecutive years only, complete predictors
model_dt <- panel[!is.na(docs_l1) & year == year_l1 + 1]
model_dt <- model_dt[is.finite(log_sjr) & is.finite(log_sjr_l1)]

model_dt <- model_dt[!is.na(nyears_l1) & !is.na(cpd_l1) & !is.na(docs_l1) & !is.na(rpd_l1)]
cat("Estimation sample:", uniqueN(model_dt$src_id), "journals,",
    nrow(model_dt), "journal-years\n")
saveRDS(model_dt, O("model_dt.rds"))
fwrite(data.table(field = all_fields, var = field_var_names), O("field_vars.csv"))
fwrite(data.table(publisher = publisher_groups, var = pub_var_names), O("publisher_vars.csv"))

## ---------- CSV exports for the Python checks (4_fe_checks.py)
m <- model_dt
fields <- fread(O("field_vars.csv")); fv <- fields$var
keep <- c("src_id","year","log_sjr","log_sjr_l1","docs_l1","cpd_l1","rpd_l1","nyears_l1","cpd","areas", fv)
fwrite(m[, ..keep], O("est_sample.csv"))
## sample including early journals (used for the density check in 4_fe_checks.py)
p <- readRDS(O("sjr_panel_full.rds"))
p <- p[src_type %in% c("journal","trade journal") & !is.na(quart) & !quart %in% c("","-")]
setorder(p, src_id, year)
p[, log_sjr := log(sjr + 0.001)]
p[, `:=`(docs_l1 = shift(docs), cpd_l1 = shift(cpd), rpd_l1 = shift(rpd), log_sjr_l1 = shift(log_sjr), year_l1 = shift(year)), by = src_id]
p <- p[!is.na(docs_l1) & !is.na(cpd_l1) & year == year_l1 + 1 & is.finite(log_sjr) & is.finite(log_sjr_l1)]
tok <- lapply(str_split(p$areas, ";"), str_trim)
for (i in seq_along(fv)) p[, (fv[i]) := as.integer(vapply(tok, function(t) fields$field[i] %in% t, logical(1)))]
p[, pre1996 := as.integer(retro_flag %in% TRUE)]
cat("full lagged sample:", nrow(p), "rows,", uniqueN(p$src_id), "journals; pre-1996 journals:", p[pre1996==1, uniqueN(src_id)], "\n")
fwrite(p[, c("src_id","year","log_sjr","log_sjr_l1","docs_l1","cpd_l1","rpd_l1","cpd","pre1996", fv), with = FALSE], O("full_sample.csv"))

## flag for the 21 large publishers (used in 4_fe_checks.py)
d <- model_dt
pv <- setdiff(grep("^pub_", names(d), value = TRUE), "pub_group")
x <- d[, .(named = as.integer(rowSums(.SD) > 0)), .SDcols = pv, by = .(src_id, year)]
fwrite(unique(x[, .(src_id, year, named)]), O("named_pub.csv"))

