## =============================================================
## 2_models.R
## Main estimates for the revised manuscript, re-based on log(CPD).
##
## Why: the slope on log(1 + CPD) is not the elasticity of SJR with
## respect to CPD. d log SJR / d log CPD = b * CPD / (1 + CPD), which
## depends on the level of CPD, so log(1 + CPD) slopes mechanically fall
## as citation rates rise and are mechanically steeper in low-citation
## fields. log(CPD) gives the elasticity directly. Journal-years with
## CPD = 0 (about 7%) are excluded from log-CPD models; log(CPD + 0.1)
## and asinh(CPD) on the full sample are reported as checks.
##
## Sections
##   A  scaling comparison: the same quantities under five scalings
##   B  M1-M4 (journal + year FE), Mundlak tests, persistence
##   C  field (M5) and publisher (M7, M7 + fields) deviations
##   D  robustness of the citation elasticity
##   E  between effects (M6), open access
##   F  density moderation (M8) with alternative clustering
##   G  sample including journals indexed from 1996 or earlier
## Run after 1_build_panel.R.
## =============================================================
suppressMessages({library(data.table); library(plm); library(lmtest); library(sandwich); library(stringr)})
source("config.R")
log_con <- file(O("2_models_log.txt"), open = "wt"); sink(log_con, split = TRUE)
set.seed(42)

d0 <- readRDS(O("model_dt.rds"))
fields <- fread(O("field_vars.csv")); pubs <- fread(O("publisher_vars.csv"))
fv <- fields$var; pv <- pubs$var
d0[, ldocs := log1p(docs_l1)]
d0[, pfield := str_trim(sapply(str_split(areas, ";"), `[`, 1))]
d <- d0[cpd_l1 > 0]
d[, lcpd := log(cpd_l1)]
setorder(d, src_id, year)
res <- list(); put <- function(k, v) res[[k]] <<- unname(v)
put("n_all", nrow(d0)); put("n_pos", nrow(d)); put("j_pos", uniqueN(d$src_id)); put("share_zero", d0[, mean(cpd_l1 == 0)])
cat("log-CPD sample:", nrow(d), "journal-years,", uniqueN(d$src_id), "journals; CPD = 0 share", round(res$share_zero, 4), "\n")

## ---------- helpers --------------------------------------------------
fe <- function(f, data = d, effect = "twoways") {
  m <- plm(f, data = data, index = c("src_id", "year"), model = "within", effect = effect)
  list(m = m, V = vcovHC(m, type = "HC1", cluster = "group"), data = data,
       r2 = unname(summary(m)$r.squared["rsq"]), n = nobs(m), J = pdim(m)$nT$n)
}
ct <- function(fit, V = fit$V) coeftest(fit$m, vcov = V)
## cluster-robust vcov for any cluster variable (rows of fit$data in model order)
vc_cl <- function(fit, cl) {
  X <- model.matrix(fit$m); u <- as.numeric(resid(fit$m))
  idx <- index(fit$m); key <- paste(idx[[1]], idx[[2]])
  dk <- paste(fit$data$src_id, fit$data$year); g <- cl[match(key, dk)]
  B <- solve(crossprod(X)); S <- rowsum(X * u, g); G <- length(unique(g)); N <- nrow(X); K <- ncol(X)
  (G / (G - 1)) * ((N - 1) / (N - K)) * B %*% crossprod(S) %*% B
}
ix <- function(x, ctrl, groups, extra = character())
  as.formula(paste("log_sjr ~", ctrl, "+", x, "+", paste(c(extra, paste0(x, ":", groups)), collapse = " + ")))
deviations <- function(fit, x, groups, labels) {
  b <- coef(fit$m); V <- fit$V; nm <- names(b); D <- fit$data
  k0 <- match(x, nm); ks <- match(paste0(x, ":", groups), nm); G <- as.matrix(D[, ..groups])
  a <- numeric(length(b)); a[k0] <- 1; a[ks] <- colMeans(G)
  out <- rbindlist(lapply(seq_along(groups), function(j) {
    g <- numeric(length(b)); g[k0] <- 1; g[ks] <- colMeans(G[G[, j] == 1, , drop = FALSE]); L <- g - a
    data.table(group = labels[j], journals = D[G[, j] == 1, uniqueN(src_id)], mean_cpd = D[G[, j] == 1, mean(cpd_l1)],
               slope = sum(g * b), dev = sum(L * b), se = sqrt(drop(t(L) %*% V %*% L)))
  }))
  out[, p := 2 * pnorm(-abs(dev / se))][, q := p.adjust(p, "BH")]
  out[, dir := fifelse(q < .05, fifelse(dev > 0, "above", "below"), "n.s.")]
  setorder(out, -dev)
  list(avg = sum(a * b), avg_se = sqrt(drop(t(a) %*% V %*% a)), tab = out)
}
year_path <- function(x, data) {
  data <- copy(data); yrs <- sort(unique(data$year))[-1]
  for (y in yrs) data[, (paste0("ly", y)) := get(x) * (year == y)]
  f <- fe(as.formula(paste("log_sjr ~ ldocs +", x, "+", paste0("ly", yrs, collapse = " + "))), data)
  b <- coef(f$m); V <- f$V
  rbindlist(lapply(c(min(data$year), yrs), function(y) {
    L <- setNames(numeric(length(b)), names(b)); L[x] <- 1; if (y != min(data$year)) L[paste0("ly", y)] <- 1
    data.table(year = y, slope = sum(L * b), se = sqrt(drop(t(L) %*% V %*% L)))
  }))
}

## ---------- A. scaling comparison -----------------------------------
scal <- list(
  list(key = "log",   lab = "log(CPD), CPD > 0", data = d,  fn = function(z) log(z)),
  list(key = "l01",   lab = "log(CPD + 0.1)",    data = d0, fn = function(z) log(z + 0.1)),
  list(key = "asinh", lab = "asinh(CPD)",        data = d0, fn = function(z) asinh(z)),
  list(key = "l1p",   lab = "log(1 + CPD)",      data = d0, fn = function(z) log1p(z)),
  list(key = "raw",   lab = "CPD (raw units)",   data = d0, fn = function(z) z))
scal_rows <- list(); paths <- list(); fdev <- list(); pdev <- list()
for (s in scal) {
  D <- copy(s$data); D[, xs := s$fn(cpd_l1)]
  m2 <- fe(log_sjr ~ ldocs + xs, D)
  pth <- year_path("xs", D); pth[, spec := s$lab]; paths[[s$key]] <- pth
  f5 <- deviations(fe(ix("xs", "ldocs", fv), D), "xs", fv, fields$field)
  f7 <- deviations(fe(ix("xs", "ldocs", pv), D), "xs", pv, pubs$publisher)
  fdev[[s$key]] <- f5$tab[, spec := s$lab]; pdev[[s$key]] <- f7$tab[, spec := s$lab]
  ## implied elasticity w.r.t. CPD at the sample's mean of dx/dlogCPD
  dxd <- switch(s$key, log = 1, l01 = D[, mean(cpd_l1 / (cpd_l1 + 0.1))], asinh = D[, mean(cpd_l1 / sqrt(1 + cpd_l1^2))],
                l1p = D[, mean(cpd_l1 / (1 + cpd_l1))], raw = D[, mean(cpd_l1)])
  scal_rows[[s$key]] <- data.table(spec = s$lab, n = m2$n, coef = coef(m2$m)["xs"], within_r2 = m2$r2,
    implied_elasticity = coef(m2$m)["xs"] * dxd,
    slope_2005 = pth[year == 2005, slope], slope_2025 = pth[year == 2025, slope],
    pct_change = pth[year == 2025, slope] / pth[year == 2005, slope] - 1,
    fields_above = f5$tab[dir == "above", .N], fields_below = f5$tab[dir == "below", .N],
    cor_field_dev_cpd = f5$tab[, cor(dev, mean_cpd)],
    pubs_above = f7$tab[dir == "above", .N], pubs_below = f7$tab[dir == "below", .N],
    ah_dir = f5$tab[group == "Arts and Humanities", dir], nursing_dir = f5$tab[group == "Nursing", dir],
    chem_dir = f5$tab[group == "Chemistry", dir])
  cat("scaling done:", s$lab, "\n")
}
scal_tab <- rbindlist(scal_rows); print(scal_tab, digits = 3)
fwrite(scal_tab, O("lc_scaling_comparison.csv")); fwrite(rbindlist(paths), O("lc_slope_paths.csv"))
fwrite(rbindlist(fdev), O("lc_field_dev_by_scaling.csv")); fwrite(rbindlist(pdev), O("lc_pub_dev_by_scaling.csv"))
## implied elasticity w.r.t. CPD by year from the log(1 + CPD) path
w <- d0[, .(share = mean(cpd_l1 / (1 + cpd_l1))), by = year]
ie <- merge(paths$l1p, w, by = "year")[, .(year, slope = slope * share, spec = "Implied elasticity from log(1 + CPD)")]
fwrite(rbind(rbindlist(paths)[, .(year, slope, spec)], ie), O("lc_slope_paths_plus.csv"))
put("ie_2005", ie[year == 2005, slope]); put("ie_2025", ie[year == 2025, slope])

## ---------- B. M1-M4, Mundlak, persistence ---------------------------
m1 <- fe(log_sjr ~ ldocs); m2 <- fe(log_sjr ~ ldocs + lcpd); m3 <- fe(log_sjr ~ ldocs + lcpd + rpd_l1)
m4 <- fe(log_sjr ~ ldocs + lcpd + log_sjr_l1)
tidy <- function(f, lab) { k <- ct(f); data.table(model = lab, term = rownames(k), est = k[, 1], se = k[, 2], p = k[, 4],
                                                  within_r2 = f$r2, n = f$n, journals = f$J) }
main <- rbind(tidy(m1, "M1"), tidy(m2, "M2"), tidy(m3, "M3"), tidy(m4, "M4")); print(main, digits = 4)
fwrite(main, O("lc_table_main.csv"))
mundlak <- function(xs, lab) {
  dd <- copy(d); for (x in xs) dd[, (paste0("mean_", x)) := mean(get(x)), by = src_id]
  l <- lm(as.formula(paste("log_sjr ~", paste(c(xs, paste0("mean_", xs)), collapse = " + "), "+ factor(year)")), dd)
  V <- vcovCL(l, cluster = ~src_id, type = "HC1"); k <- paste0("mean_", xs); b <- coef(l)[k]
  W <- drop(t(b) %*% solve(V[k, k]) %*% b); data.table(model = lab, chi2 = W, df = length(k), p = pchisq(W, length(k), lower.tail = FALSE))
}
haus <- rbind(mundlak("ldocs", "M1"), mundlak(c("ldocs", "lcpd"), "M2"), mundlak(c("ldocs", "lcpd", "rpd_l1"), "M3"))
print(haus); fwrite(haus, O("lc_table_mundlak.csv"))
f4 <- log_sjr ~ ldocs + lcpd + log_sjr_l1
rho_fe <- coef(m4$m)["log_sjr_l1"]; rho_pool <- coef(lm(update(f4, . ~ . + factor(year)), d))["log_sjr_l1"]
d[, `:=`(pos = seq_len(.N) - 1L, Tn = .N), by = src_id]
## split-panel jackknife of both the persistence and the CPD coefficient, so that the
## long-run elasticity b / (1 - rho) uses bias-corrected values of both
jk <- function(D) {
  cf <- function(H) coef(plm(f4, H, index = c("src_id", "year"), model = "within", effect = "twoways"))[c("log_sjr_l1", "lcpd")]
  r <- cf(D); h <- sapply(list(D[pos < Tn %/% 2], D[pos >= Tn %/% 2]), cf)
  j <- 2 * r - rowMeans(h)
  c(full = unname(r[1]), h1 = unname(h[1, 1]), h2 = unname(h[1, 2]), jk = unname(j[1]),
    b_jk = unname(j[2]), lr_jk = unname(j[2] / (1 - j[1])))
}
j0 <- jk(d); ids <- unique(d$src_id)
boot <- replicate(100, { s <- data.table(src_id = sample(ids, replace = TRUE))[, bid := .I]
  b <- merge(s, d, by = "src_id", allow.cartesian = TRUE); b[, src_id := bid]; jk(b)[c("jk", "lr_jk")] })
put("rho_fe", rho_fe); put("rho_pool", rho_pool); put("rho_h1", j0["h1"]); put("rho_h2", j0["h2"]); put("rho_jk", j0["jk"]); put("rho_jk_se", sd(boot["jk", ]))
put("b_jk", j0["b_jk"]); put("lr_jk", j0["lr_jk"]); put("lr_jk_se", sd(boot["lr_jk", ]))
put("T_mean", d[, .N, by = src_id][, mean(N)])
cat(sprintf("persistence: pooled %.3f FE %.3f halves %.3f %.3f jackknife %.3f (SE %.3f)\n", rho_pool, rho_fe, j0["h1"], j0["h2"], j0["jk"], sd(boot["jk", ])))

## ---------- C. field and publisher deviations (log CPD) ---------------
m5 <- fe(ix("lcpd", "ldocs", fv)); dv5 <- deviations(m5, "lcpd", fv, fields$field)
m7 <- fe(ix("lcpd", "ldocs", pv)); dv7 <- deviations(m7, "lcpd", pv, pubs$publisher)
m7f <- fe(ix("lcpd", "ldocs", c(pv, fv))); dv7f <- deviations(m7f, "lcpd", pv, pubs$publisher)
for (k in c("5", "7", "7f")) { dv <- get(paste0("dv", k)); fwrite(dv$tab, O(paste0("lc_dev_m", k, ".csv")))
  put(paste0("avg", k), dv$avg); put(paste0("avg", k, "_se"), dv$avg_se)
  put(paste0("above", k), dv$tab[dir == "above", .N]); put(paste0("below", k), dv$tab[dir == "below", .N]) }
put("cor5", dv5$tab[, cor(dev, mean_cpd)]); put("cor7", dv7$tab[, cor(dev, mean_cpd)])
put("ref7_j", d[rowSums(d[, ..pv]) == 0, uniqueN(src_id)]); put("ref7_cpd", d[rowSums(d[, ..pv]) == 0, mean(cpd_l1)])
put("in7_cpd", d[rowSums(d[, ..pv]) > 0, mean(cpd_l1)]); put("mean_fields", mean(rowSums(as.matrix(d[, ..fv]))))
print(dv5$tab, digits = 3); print(dv7$tab, digits = 3); print(dv7f$tab, digits = 3)

## ---------- D. robustness of the elasticity -----------------------------
rob <- list(); addr <- function(lab, f, term) { k <- ct(f); rob[[lab]] <<- data.table(spec = lab, est = k[term, 1], se = k[term, 2], p = k[term, 4], n = f$n) }
addr("M2: journal + year FE, log CPD", m2, "lcpd")
addr("Journal FE only, NYears added", fe(log_sjr ~ ldocs + nyears_l1 + lcpd, effect = "individual"), "lcpd")
D01 <- copy(d0)[, x := log(cpd_l1 + 0.1)]; addr("log(CPD + 0.1), all journal-years", fe(log_sjr ~ ldocs + x, D01), "x")
DA <- copy(d0)[, x := asinh(cpd_l1)]; addr("asinh(CPD), all journal-years", fe(log_sjr ~ ldocs + x, DA), "x")
h3 <- d0[, .(src_id, year = year + 3L, cpd_l3 = cpd)]
d3 <- merge(d, h3, by = c("src_id", "year"))[cpd_l3 > 0][, lcpd3 := log(cpd_l3)]
addr("Lag-3 sample: CPD at t-1", fe(log_sjr ~ ldocs + lcpd, d3), "lcpd")
addr("Lag-3 sample: CPD at t-3 (no window overlap)", fe(log_sjr ~ ldocs + lcpd3, d3), "lcpd3")
addr("Outcome years 2005-2014", fe(log_sjr ~ ldocs + lcpd, d[year <= 2014]), "lcpd")
addr("Outcome years 2015-2025", fe(log_sjr ~ ldocs + lcpd, d[year >= 2015]), "lcpd")
rob_tab <- rbindlist(rob); print(rob_tab, digits = 4); fwrite(rob_tab, O("lc_table_robustness.csv"))
d[, late := as.integer(year >= 2015)]; lt <- ct(fe(log_sjr ~ ldocs + lcpd + lcpd:late))
put("late_int", lt["lcpd:late", 1]); put("late_int_se", lt["lcpd:late", 2]); put("late_int_p", lt["lcpd:late", 4])

## ---------- E. between effects (OA) -------------------------------------
be <- d[, c(list(log_sjr = mean(log_sjr), ldocs = mean(ldocs), nyears = mean(nyears_l1), lcpd = mean(lcpd),
                 oa = first(oa_status)), lapply(.SD, max)), by = src_id, .SDcols = fv]
## the scaling variants use the same journals and the same controls; only the citation term changes
be[, `:=`(l1p = d[, mean(log1p(cpd_l1)), by = src_id]$V1, cpd = d[, mean(cpd_l1), by = src_id]$V1)]
fld <- paste(fv, collapse = " + ")
bfit <- function(rhs, data, lab) { l <- lm(as.formula(paste("log_sjr ~", rhs)), data); k <- coeftest(l, vcov = vcovHC(l, type = "HC3"))
  list(lab = lab, k = k, r2 = summary(l)$r.squared, n = nobs(l)) }
bm <- bfit(paste("ldocs + nyears + lcpd + oa +", fld), be, "M6: mean log CPD + subject areas")
bl <- list(bm, bfit("ldocs + nyears + lcpd + oa", be, "Mean log CPD, no subject areas"),
           bfit(paste("ldocs + nyears + oa +", fld), be, "No CPD, subject areas"),
           bfit(paste("ldocs + nyears + l1p + oa +", fld), be, "Mean log(1 + CPD), subject areas"),
           bfit(paste("ldocs + nyears + cpd + oa +", fld), be, "Mean raw CPD, subject areas"))
be_tab <- rbindlist(lapply(bl, function(b) data.table(spec = b$lab, oa = b$k["oa", 1], oa_se = b$k["oa", 2], oa_p = b$k["oa", 4], r2 = b$r2, n = b$n)))
print(be_tab, digits = 4); fwrite(be_tab, O("lc_table_between_oa.csv"))
fwrite(data.table(term = rownames(bm$k), bm$k[, c(1, 2, 4)]), O("lc_table_between_main.csv"))
put("be_n", bm$n); put("be_r2", bm$r2)

## ---------- F. density moderation (M8) ------------------------------------
G <- as.matrix(d[, ..fv]); nf <- rowSums(G); d[, lrpd := log(pmax(rpd_l1, 0.5))]
fymean <- function(col) { M <- sapply(seq_along(fv), function(j) { fy <- d[G[, j] == 1, .(v = mean(get(col))), by = year]
  fy$v[match(d$year, fy$year)] }); M[is.na(M)] <- 0; rowSums(G * M) / nf }
d[, dens := fymean("lcpd")]; d[, refd := fymean("lrpd")]
d[, `:=`(dens_c = dens - mean(dens), refd_c = refd - mean(refd), lsq = (lcpd - mean(lcpd))^2)]
put("dens_sd", sd(d$dens)); put("refd_sd", sd(d$refd))
yrs <- sort(unique(d$year))[-1]; for (y in yrs) d[, (paste0("ly", y)) := lcpd * (year == y)]; ly <- paste0("ly", yrs)
d[, fyc := paste(pfield, year)]
m8rows <- list()
for (dv in c("dens_c", "refd_c")) for (spec in c("pooled", "field", "strict")) {
  extra <- switch(spec, pooled = character(), field = paste0("lcpd:", fv), strict = c(paste0("lcpd:", fv), ly, "lsq"))
  f <- fe(as.formula(paste("log_sjr ~ ldocs + lcpd +", dv, "+", paste0("lcpd:", dv), if (length(extra)) paste("+", paste(extra, collapse = " + ")))))
  term <- paste0("lcpd:", dv); b <- coef(f$m)[term]
  se_j <- sqrt(f$V[term, term])
  if (spec == "strict") {
    Vf <- vc_cl(f, d$pfield); Vy <- vc_cl(f, as.character(d$year)); Vjy <- vc_cl(f, paste(d$src_id, d$year)); Vj <- vc_cl(f, as.character(d$src_id))
    Vfy <- vc_cl(f, d$fyc); V2 <- Vj + Vy - Vjy
    se_f <- sqrt(Vf[term, term]); se_2 <- sqrt(max(V2[term, term], 0)); se_fy <- sqrt(Vfy[term, term])
  } else { se_f <- se_2 <- se_fy <- NA_real_ }
  m8rows[[paste(dv, spec)]] <- data.table(density = ifelse(dv == "dens_c", "Citation density (mean log CPD)", "Referencing density (mean log RPD)"),
    spec = spec, est = b, se_journal = se_j, se_field = se_f, se_twoway = se_2, se_field_year = se_fy, within_r2 = f$r2)
  cat("M8", dv, spec, round(b, 4), round(se_j, 4), round(se_f, 4), "\n")
}
m8 <- rbindlist(m8rows); print(m8, digits = 4); fwrite(m8, O("lc_table_density.csv"))

## ---------- G. including journals indexed from 1996 or earlier ------------
p <- readRDS(O("sjr_panel_full.rds"))
p <- p[src_type %in% c("journal", "trade journal") & !is.na(quart) & !quart %in% c("", "-")]
setorder(p, src_id, year); p[, log_sjr := log(sjr + 0.001)]
p[, `:=`(docs_l1 = shift(docs), cpd_l1 = shift(cpd), year_l1 = shift(year)), by = src_id]
p <- p[!is.na(docs_l1) & !is.na(cpd_l1) & year == year_l1 + 1 & is.finite(log_sjr) & cpd_l1 > 0]
p[, `:=`(ldocs = log1p(docs_l1), lcpd = log(cpd_l1))]
fg <- fe(log_sjr ~ ldocs + lcpd, p); pg <- year_path("lcpd", p)
put("full_n", fg$n); put("full_j", fg$J); put("full_m2", coef(fg$m)["lcpd"]); put("full_m2_se", sqrt(fg$V["lcpd", "lcpd"]))
put("full_2005", pg[year == 2005, slope]); put("full_2025", pg[year == 2025, slope])
put("full_pre1996_j", p[retro_flag %in% TRUE, uniqueN(src_id)])
cat("including pre-1996:", fg$n, "obs", fg$J, "journals; elasticity", round(coef(fg$m)["lcpd"], 3), "; 2005", round(pg[year == 2005, slope], 3), "2025", round(pg[year == 2025, slope], 3), "\n")

jsonlite::write_json(res, O("lc_numbers.json"), auto_unbox = TRUE, digits = 10)
str(res); sink()
