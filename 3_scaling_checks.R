## =============================================================
## 3_scaling_checks.R
## Scale invariance of field and publisher conclusions (Table 6), the
## large-publisher gap by scaling (Table 7), conversion of coefficients
## into elasticities (Figure 3), descriptives (Table 2), and the
## same-sample and standard-error checks. Run after 2_models.R.
## =============================================================
suppressMessages({library(data.table); library(plm); library(lmtest); library(sandwich); library(stringr)})
source("config.R")

## ---------- scale invariance (Table 6)
## =============================================================
## 11_scale_invariance.R
## For each subject area and publisher: direction of its deviation
## from the average slope (BH q < .05) under each of the five scalings,
## and a summary of which conclusions hold under all of them.
## Run after 09_logcpd_models.R.
## =============================================================
logs <- c("log(CPD), CPD > 0", "log(CPD + 0.1)", "asinh(CPD)", "log(1 + CPD)")
for (k in c("field", "pub")) {
  d <- fread(O(sprintf("lc_%s_dev_by_scaling.csv", k)))
  w <- dcast(d, group ~ spec, value.var = "dir")
  allv <- names(w)[-1]
  w[, all5 := apply(.SD, 1, function(r) if (length(unique(r)) == 1 && r[1] != "n.s.") r[1] else ""), .SDcols = allv]
  w[, all4log := apply(.SD, 1, function(r) if (length(unique(r)) == 1 && r[1] != "n.s.") r[1] else ""), .SDcols = logs]
  w[, flips := apply(.SD, 1, function(r) "above" %in% r && "below" %in% r), .SDcols = logs]
  fwrite(w, O(sprintf("lc_invariance_%s.csv", k)))
  cat(k, ": all five", w[all5 != "", .N], "| all four log-type", w[all4log != "", .N], "| opposite signs across log-type", w[flips == TRUE, .N], "\n")
}

## ---------- large-publisher gap by scaling (Table 7)
suppressMessages({library(data.table); library(plm); library(sandwich); library(lmtest)})
d0 <- readRDS(O("model_dt.rds")); d0[, ldocs := log1p(docs_l1)]
pv <- fread(O("publisher_vars.csv"))$var; fv <- fread(O("field_vars.csv"))$var
d0[, big := as.integer(rowSums(.SD) > 0), .SDcols = pv]
sc <- list(log = list(d0[cpd_l1 > 0], function(z) log(z)), l01 = list(d0, function(z) log(z + 0.1)),
           asinh = list(d0, asinh), l1p = list(d0, log1p), raw = list(d0, identity))
out <- rbindlist(lapply(names(sc), function(k) {
  D <- copy(sc[[k]][[1]]); D[, x := sc[[k]][[2]](cpd_l1)]
  r <- lapply(c(FALSE, TRUE), function(withf) {
    f <- as.formula(paste("log_sjr ~ ldocs + x + x:big", if (withf) paste("+", paste0("x:", fv, collapse = " + ")) else ""))
    m <- plm(f, D, index = c("src_id", "year"), model = "within", effect = "twoways")
    kk <- coeftest(m, vcov = vcovHC(m, type = "HC1", cluster = "group"))
    data.table(scaling = k, fields = withf, ref = kk["x", 1], diff = kk["x:big", 1], se = kk["x:big", 2], ratio = kk["x:big", 1] / kk["x", 1])
  }); rbindlist(r)
}))
print(out, digits = 3); fwrite(out, O("lc_bigpub_by_scaling.csv"))

## ---------- large-publisher gap in log(CPD) (added to lc_numbers.json)
d <- readRDS(O("model_dt.rds"))[cpd_l1 > 0]; d[, `:=`(ldocs = log1p(docs_l1), lcpd = log(cpd_l1))]
pv <- fread(O("publisher_vars.csv"))$var
d[, big := as.integer(rowSums(.SD) > 0), .SDcols = pv]
m <- plm(log_sjr ~ ldocs + lcpd + lcpd:big, d, index = c("src_id", "year"), model = "within", effect = "twoways")
k <- coeftest(m, vcov = vcovHC(m, type = "HC1", cluster = "group")); print(k)
x <- jsonlite::fromJSON(O("lc_numbers.json"))
x$ref_slope <- k["lcpd", 1]; x$big_diff <- k["lcpd:big", 1]; x$big_diff_se <- k["lcpd:big", 2]; x$ref_obs_share <- d[, mean(big == 0)]
jsonlite::write_json(x, O("lc_numbers.json"), auto_unbox = TRUE, digits = 10)
cat("ref obs share", d[, mean(big == 0)], "\n")

## ---------- converted elasticities (Figure 3)
## 12_converted_elasticity.R  (run after 09_logcpd_models.R and 11_scale_invariance.R)
## Convert each group's log(1+CPD), asinh and raw slopes to elasticities at the group's own CPD,
## and compare with the log(CPD) slopes. Tests whether scale dependence is the transformation's derivative.
d0 <- readRDS(O("model_dt.rds")); fv <- grep("^field_", names(d0), value = TRUE)
pv <- setdiff(grep("^pub_", names(d0), value = TRUE), "pub_group")
fd <- fread(O("lc_field_dev_by_scaling.csv")); pd <- fread(O("lc_pub_dev_by_scaling.csv"))
fac <- function(rows) list(l1p = mean(rows$cpd_l1 / (1 + rows$cpd_l1)), asinh = mean(rows$cpd_l1 / sqrt(1 + rows$cpd_l1^2)),
                           l01 = mean(rows$cpd_l1 / (rows$cpd_l1 + 0.1)), raw = mean(rows$cpd_l1))
key <- c(l1p = "log(1 + CPD)", asinh = "asinh(CPD)", l01 = "log(CPD + 0.1)", raw = "CPD (raw units)")
conv <- function(tab, vars, labs) {
  out <- rbindlist(lapply(seq_along(vars), function(k) { r <- d0[get(vars[k]) == 1]; fc <- fac(r)
    rbindlist(lapply(names(key), function(s) data.table(group = labs[k], sk = s, factor = fc[[s]]))) }))
  m <- merge(out, tab[, .(group, spec, slope)], by.x = c("group"), by.y = "group", allow.cartesian = TRUE)[spec == key[sk]]
  m[, el := slope * factor]
  lg <- tab[spec == "log(CPD), CPD > 0", .(group, el_log = slope)]
  m <- merge(m, lg, by = "group")
  m[, `:=`(dev_c = el - mean(el), dev_log = el_log - mean(el_log), sdev = slope - mean(slope)), by = sk]
  list(sum = m[, .(r_raw_slope = cor(slope, el_log), r_converted = cor(el, el_log), mean_abs_gap = mean(abs(el - el_log)),
                   sign_agree_raw = mean(sign(sdev) == sign(dev_log)), sign_agree_conv = mean(sign(dev_c) == sign(dev_log)), n = .N), by = sk],
       detail = m)
}
flab <- sub("^field_", "", fv); flab <- fd$group[match(gsub("[^A-Za-z]", "", flab), gsub("[^A-Za-z]", "", fd$group))]
plab <- pd$group[match(gsub("[^A-Za-z]", "", sub("^pub_", "", pv)), gsub("[^A-Za-z]", "", pd$group))]
stopifnot(!anyNA(flab), !anyNA(plab))
rf <- conv(fd, fv, flab); rp <- conv(pd, pv, plab)
print(rf$sum); print(rp$sum)
fwrite(rbind(rf$sum[, level := "field"], rp$sum[, level := "publisher"]), O("lc_converted_elasticity.csv"))
fwrite(rbind(rf$detail[, level := "field"], rp$detail[, level := "publisher"]), O("lc_converted_elasticity_detail.csv"))
## derivative factors over time and for the named-publisher vs residual groups
yr <- d0[year %in% c(2005, 2025), .(mean_cpd = mean(cpd_l1), fac_l1p = mean(cpd_l1 / (1 + cpd_l1))), by = year]
big <- d0[, .(fac_l1p = mean(cpd_l1 / (1 + cpd_l1)), fac_asinh = mean(cpd_l1 / sqrt(1 + cpd_l1^2))), by = .(named = rowSums(d0[, ..pv]) > 0)]
print(yr); print(big)
bp <- fread(O("lc_bigpub_by_scaling.csv"))[fields == FALSE]
cvb <- function(s, fcol) { r <- bp[scaling == s]; fr <- big[named == FALSE][[fcol]]; fb <- big[named == TRUE][[fcol]]
  c(ref = r$ref * fr, named = (r$ref + r$diff) * fb) }
cb <- rbind(l1p = cvb("l1p", "fac_l1p"), asinh = cvb("asinh", "fac_asinh")); print(cb)
jsonlite::write_json(list(cpd_2005 = yr[year == 2005, mean_cpd], cpd_2025 = yr[year == 2025, mean_cpd],
  fac_2005 = yr[year == 2005, fac_l1p], fac_2025 = yr[year == 2025, fac_l1p],
  big_conv_l1p_ref = cb["l1p", "ref"], big_conv_l1p_named = cb["l1p", "named"],
  big_conv_asinh_ref = cb["asinh", "ref"], big_conv_asinh_named = cb["asinh", "named"]), O("lc_conv_numbers.json"), auto_unbox = TRUE, digits = 10)

## ---------- descriptives (Table 2)
## 13_descriptives_logcpd.R  (run after 09_logcpd_models.R and 11_scale_invariance.R)
d0 <- readRDS(O("model_dt.rds")); d <- d0[cpd_l1 > 0]

v <- list(SJR = exp(d$log_sjr) - 0.001, Documents = d$docs_l1, CPD = d$cpd_l1, RPD = d$rpd_l1, NYears = d$nyears_l1)
out <- rbindlist(lapply(names(v), function(n) { x <- v[[n]]; data.table(variable = n, mean = mean(x), sd = sd(x),
  p25 = quantile(x, .25), median = median(x), p75 = quantile(x, .75), max = max(x)) }))
print(out); fwrite(out, O("lc_table_descriptives.csv"))
setorder(d, src_id, year)
inc <- d[, .(dd = diff(nyears_l1), dy = diff(year)), by = src_id][dy == 1]
oa_j <- d[, .(oa = oa_status[1]), by = src_id][, mean(oa)]
ly <- d[, .(m = mean(log_sjr)), by = year][order(year)]
# zero-CPD journals: describe
z <- d0[cpd_l1 == 0]; zf <- d0[, .(any0 = any(cpd_l1 == 0), all0 = all(cpd_l1 == 0)), by = src_id]
res <- list(nyears_inc1 = inc[, mean(dd == 1)], oa_share_j = oa_j, logsjr_2005 = ly[year == 2005, m], logsjr_2025 = ly[year == 2025, m],
  zero_j_any = zf[, sum(any0)], zero_j_all = zf[, sum(all0)], zero_mean_sjr = z[, mean(exp(log_sjr) - .001)], pos_mean_sjr = d[, mean(exp(log_sjr) - .001)],
  zero_share_ah = NA)
## zero shares by period and field (written below)
print(res); jsonlite::write_json(res, O("lc_desc_numbers.json"), auto_unbox = TRUE, digits = 10)
fv <- grep("^field_", names(d0), value = TRUE)
res$zero_early <- d0[year <= 2014, mean(cpd_l1 == 0)]; res$zero_late <- d0[year >= 2015, mean(cpd_l1 == 0)]
res$zero_ah <- d0[field_Arts.and.Humanities == 1, mean(cpd_l1 == 0)]
jsonlite::write_json(res, O("lc_desc_numbers.json"), auto_unbox = TRUE, digits = 10)

## ---------- same-sample comparison and converted elasticities with SEs
## =============================================================
## 16_review_checks.R
## Responses to review comments:
##  (1) elasticity of SJR itself vs. of SJR + 0.001 (outcome constant)
##  (2) all five scalings on the SAME positive-CPD sample, so that scaling
##      is separated from sample composition (zeros)
##  (3) converted elasticities with standard errors, computed by averaging
##      observation-level elasticities (handles multi-field journals),
##      and confidence intervals for the agreement correlations
## Run after 09_logcpd_models.R.
## =============================================================
src <- readLines("2_models.R")
cut <- grep("^## ---------- A. scaling comparison", src)[1] - 1
src <- src[seq_len(cut)]; src <- src[!grepl("sink\\(|log_con", src)]
eval(parse(text = src))
out <- list(); put2 <- function(k, v) out[[k]] <<- unname(v)

## (1) outcome constant
put2("sjr_min", d[, min(exp(log_sjr) - 0.001)])
put2("const_factor_mean", d[, mean(exp(log_sjr) / (exp(log_sjr) - 0.001))])
d[, log_sjr_exact := log(exp(log_sjr) - 0.001)]
m_exact <- fe(log_sjr_exact ~ ldocs + lcpd, d)
put2("m2_exact", coef(m_exact$m)["lcpd"]); put2("m2_exact_se", sqrt(m_exact$V["lcpd", "lcpd"]))

## observation-level conversion of group slopes: elasticity_i = slope_i * factor_i
conv_groups <- function(fit, x, groups, labels, fac) {
  b <- coef(fit$m); V <- fit$V; nm <- names(b); D <- fit$data
  k0 <- match(x, nm); ks <- match(paste0(x, ":", groups), nm); G <- as.matrix(D[, ..groups])
  rbindlist(lapply(seq_along(groups), function(j) {
    rows <- G[, j] == 1
    ## average of factor-weighted membership vectors: E[factor_i * (1, G_i)]
    w <- fac[rows]
    g <- numeric(length(b)); g[k0] <- mean(w); g[ks] <- colMeans(G[rows, , drop = FALSE] * w)
    data.table(group = labels[j], el = sum(g * b), el_se = sqrt(drop(t(g) %*% V %*% g)))
  }))
}
## sample-average slope (fac = 1) or elasticity: the reference of Equation (3)
conv_all <- function(fit, x, groups, fac) {
  b <- coef(fit$m); V <- fit$V; nm <- names(b); D <- fit$data
  k0 <- match(x, nm); ks <- match(paste0(x, ":", groups), nm); G <- as.matrix(D[, ..groups])
  g <- numeric(length(b)); g[k0] <- mean(fac); g[ks] <- colMeans(G * fac)
  c(est = sum(g * b), se = sqrt(drop(t(g) %*% V %*% g)))
}
## group slopes and elasticities as deviations from the sample average, before and after conversion
group_devs <- function(fit, x, groups, labels, fac, key, lev) {
  one <- rep(1, length(fac))
  s <- conv_groups(fit, x, groups, labels, one); e <- conv_groups(fit, x, groups, labels, fac)
  data.table(group = labels, level = lev, sk = key, slope = s$el, slope_dev = s$el - conv_all(fit, x, groups, one)["est"],
             el = e$el, el_se = e$el_se, el_dev = e$el - conv_all(fit, x, groups, fac)["est"])
}
fisher <- function(r, n) { z <- atanh(r); s <- 1 / sqrt(n - 3); tanh(c(z - 1.96 * s, z + 1.96 * s)) }

scal <- list(
  list(key = "log",   lab = "log(CPD)",       fn = function(z) log(z),          fac = function(z) rep(1, length(z))),
  list(key = "l01",   lab = "log(CPD + 0.1)", fn = function(z) log(z + 0.1),    fac = function(z) z / (z + 0.1)),
  list(key = "asinh", lab = "asinh(CPD)",     fn = function(z) asinh(z),        fac = function(z) z / sqrt(1 + z^2)),
  list(key = "l1p",   lab = "log(1 + CPD)",   fn = function(z) log1p(z),        fac = function(z) z / (1 + z)),
  list(key = "raw",   lab = "CPD (raw units)",fn = function(z) z,               fac = function(z) z))
rows <- list(); fconv <- list(); pconv <- list(); gdev <- list()
for (s in scal) {
  D <- copy(d); D[, xs := s$fn(cpd_l1)]; fac <- s$fac(D$cpd_l1)
  m2 <- fe(log_sjr ~ ldocs + xs, D)
  pth <- year_path("xs", D)
  f5 <- fe(ix("xs", "ldocs", fv), D); f7 <- fe(ix("xs", "ldocs", pv), D)
  dv5 <- deviations(f5, "xs", fv, fields$field)$tab
  fconv[[s$key]] <- conv_groups(f5, "xs", fv, fields$field, fac)[, sk := s$key]
  pconv[[s$key]] <- conv_groups(f7, "xs", pv, pubs$publisher, fac)[, sk := s$key]
  gdev[[paste0(s$key, "f")]] <- group_devs(f5, "xs", fv, fields$field, fac, s$key, "field")
  gdev[[paste0(s$key, "p")]] <- group_devs(f7, "xs", pv, pubs$publisher, fac, s$key, "publisher")
  ## conversion of the year path at each year's mean factor
  fy <- D[, .(f = mean(s$fac(cpd_l1))), by = year]
  pth <- merge(pth, fy, by = "year")
  rows[[s$key]] <- data.table(spec = s$lab, n = m2$n, coef = coef(m2$m)["xs"],
    elasticity = coef(m2$m)["xs"] * mean(fac),
    slope_2005 = pth[year == 2005, slope], slope_2025 = pth[year == 2025, slope],
    coef_change = pth[year == 2025, slope] / pth[year == 2005, slope] - 1,
    el_2005 = pth[year == 2005, slope * f], el_2025 = pth[year == 2025, slope * f],
    el_change = pth[year == 2025, slope * f] / pth[year == 2005, slope * f] - 1,
    fields_above = dv5[dir == "above", .N], fields_below = dv5[dir == "below", .N],
    cor_field_dev_cpd = dv5[, cor(dev, mean_cpd)])
  cat("same-sample scaling done:", s$lab, "\n")
}
tab <- rbindlist(rows); print(tab, digits = 3)
fwrite(tab, O("lc_same_sample_scaling.csv"))
FC <- rbindlist(fconv); PC <- rbindlist(pconv)
fwrite(rbind(FC[, level := "field"], PC[, level := "publisher"]), O("lc_converted_with_se.csv"))
agree <- function(C, n) {
  base <- C[sk == "log", .(group, el_log = el)]
  C[sk != "log"][base, on = "group"][, {
    r <- cor(el, el_log); ci <- fisher(r, .N)
    .(r = r, r_lo = ci[1], r_hi = ci[2], mean_abs_gap = mean(abs(el - el_log)),
      mean_se = mean(el_se), sign_agree = mean(sign(el - mean(el)) == sign(el_log - mean(el_log))))
  }, by = sk]
}
AF <- agree(FC); AP <- agree(PC); print(AF); print(AP)
fwrite(rbind(AF[, level := "field"], AP[, level := "publisher"]), O("lc_same_sample_agreement.csv"))

## main comparison (Figure 3B-C): deviations from the sample average, before and after conversion
GD <- rbindlist(gdev); fwrite(GD, O("lc_group_devs.csv"))
base <- GD[sk == "log", .(group, level, dev_log = el_dev, el_log = el)]
GA <- GD[sk != "log"][base, on = c("group", "level")][, {
  rb <- cor(slope_dev, dev_log); ra <- cor(el_dev, dev_log); ci <- fisher(ra, .N)
  .(r_before = rb, r_after = ra, r_after_lo = ci[1], r_after_hi = ci[2],
    sign_before = sum(sign(slope_dev) == sign(dev_log)), sign_after = sum(sign(el_dev) == sign(dev_log)), n = .N,
    level_gap = mean(abs(el - el_log)), mean_se = mean(el_se))
}, by = .(level, sk)]
print(GA); fwrite(GA, O("lc_group_agreement.csv"))
jsonlite::write_json(out, O("lc_review_numbers.json"), auto_unbox = TRUE, digits = 10)
print(out)
