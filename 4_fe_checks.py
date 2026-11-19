"""4_fe_checks.py
Three checks estimated with pyfixest. Run after 1_build_panel.R and 2_models.R.
  (a) primary field x year fixed effects (Table 8, Section 3.5)
  (b) composition checks: year slopes on journals observed in every year; field deviations with early journals
  (c) journal size: elasticity by size quartile, document weights, larger journals, large-publisher gap net of size
"""
import os, json, gc, warnings
import numpy as np, pandas as pd, pyfixest as pf
warnings.filterwarnings("ignore")
ROOT = os.path.join(os.environ.get("DATA_DIR", "data"), "built") + os.sep

# ---------- (a) primary field x year fixed effects
d = pd.read_csv(ROOT + "est_sample.csv", usecols=lambda c: c in ("src_id", "year", "log_sjr", "docs_l1", "cpd_l1", "areas"))
d = d.merge(pd.read_csv(ROOT + "named_pub.csv"), on=["src_id", "year"], how="left")
d["named"] = d.named.fillna(0)
d["ldocs"] = np.log1p(d.docs_l1); d["l1p"] = np.log1p(d.cpd_l1)
d["pfield"] = d.areas.astype(str).str.split(";").str[0].str.strip()
d["fe_id"] = d.src_id.astype(str)
pos = d[d.cpd_l1 > 0].copy(); pos["lcpd"] = np.log(pos.cpd_l1)
out = {}
def est(fml, data, term, key):
    m = pf.feols(fml, data=data, vcov={"CRV1": "fe_id"}); t = m.tidy().loc[term]
    out[key] = {"est": float(t.Estimate), "se": float(t["Std. Error"]), "n": int(m._N)}
    print(key, out[key], flush=True); return m
est("log_sjr ~ ldocs + lcpd | fe_id + year", pos, "lcpd", "m2_base")
est("log_sjr ~ ldocs + lcpd | fe_id + pfield^year", pos, "lcpd", "m2_fy")
m = pf.feols("log_sjr ~ ldocs + lcpd | fe_id + pfield^year", data=pos, vcov={"CRV1": "pfield"})
out["m2_fy"]["se_field"] = float(m.tidy().loc["lcpd", "Std. Error"])
# year-specific slopes with field-year FE
yrs = sorted(pos.year.unique())[1:]
for y in yrs: pos[f"ly{y}"] = pos.lcpd * (pos.year == y)
m = pf.feols("log_sjr ~ ldocs + lcpd + " + " + ".join(f"ly{y}" for y in yrs) + " | fe_id + pfield^year", data=pos, vcov={"CRV1": "fe_id"})
b = m.coef(); out["fy_2005"] = float(b["lcpd"]); out["fy_2025"] = float(b["lcpd"] + b["ly2025"])
# log(1 + CPD) on all journal-years: does the apparent decline survive field-year FE?
for y in yrs: d[f"my{y}"] = d.l1p * (d.year == y)
m = pf.feols("log_sjr ~ ldocs + l1p + " + " + ".join(f"my{y}" for y in yrs) + " | fe_id + pfield^year", data=d, vcov={"CRV1": "fe_id"})
b = m.coef(); out["fy_l1p_2005"] = float(b["l1p"]); out["fy_l1p_2025"] = float(b["l1p"] + b["my2025"])
# large-publisher gap
pos["lcpd_named"] = pos.lcpd * pos.named
est("log_sjr ~ ldocs + lcpd + lcpd_named | fe_id + pfield^year", pos, "lcpd_named", "big_fy")
print({k: out[k] for k in ("fy_2005", "fy_2025", "fy_l1p_2005", "fy_l1p_2025")})
out["n_pfield"] = int(pos.pfield.nunique())
json.dump(out, open(ROOT + "lc_field_year_fe.json", "w"), indent=1)

# ---------- (b) composition checks (Section 3.5)
def prep(path):
    d = pd.read_csv(path); fv = [c for c in d.columns if c.startswith("field_")]
    d["fe_id"] = d.src_id.astype(str); d["ldocs"] = np.log1p(d.docs_l1)
    d = d[d.cpd_l1 > 0].copy(); d["lcpd"] = np.log(d.cpd_l1)
    return d, fv
full, fvf = prep(ROOT + "full_sample.csv")
out = {}
def year_slopes(d):
    yrs = sorted(d.year.unique())[1:]
    for y in yrs: d[f"ly{y}"] = d.lcpd * (d.year == y)
    b = pf.feols("log_sjr ~ ldocs + lcpd + " + " + ".join(f"ly{y}" for y in yrs) + " | fe_id + year", data=d, vcov={"CRV1": "fe_id"}).coef()
    return float(b["lcpd"]), float(b["lcpd"] + b["ly2025"])
nyr = pos.groupby("src_id").year.nunique()
bal = pos[pos.src_id.isin(nyr[nyr == pos.year.nunique()].index)].copy()
out["balanced_j"] = int(bal.src_id.nunique()); out["balanced_2005"], out["balanced_2025"] = year_slopes(bal)
def field_devs(d, fv):
    m = pf.feols("log_sjr ~ ldocs + lcpd + " + " + ".join(f"lcpd:{f}" for f in fv) + " | fe_id + year", data=d, vcov={"CRV1": "fe_id"})
    b = np.array([m.coef().get(f"lcpd:{f}", 0.0) for f in fv]); G = d[fv].values.astype(float); Gbar = G.mean(0)
    return np.array([b @ (G[G[:, j] == 1].mean(0) - Gbar) for j in range(len(fv))])   # Equation (3)
main = full[full.src_id.isin(set(pos.src_id))]
early = full[full.src_id.isin(set(pos.src_id)) | (full.pre1996 == 1)]
a, e = field_devs(main.copy(), fvf), field_devs(early.copy(), fvf)
out["early_r"] = float(np.corrcoef(a, e)[0, 1]); out["early_sign_agree"] = int((np.sign(a) == np.sign(e)).sum())
print(out); json.dump(out, open(ROOT + "lc_composition.json", "w"), indent=1)

# ---------- (c) journal size (Table 8, Section 3.5)
out = {}
pos["jsize"] = pos.groupby("src_id").docs_l1.transform("median")
pos["ls"] = np.log1p(pos.jsize); pos["ls"] -= pos.ls.mean(); pos["lcpd_ls"] = pos.lcpd * pos.ls
pos["sq"] = pd.qcut(pos.jsize, 4, labels=False)
for q in range(4):
    s = pos[pos.sq == q]; est("log_sjr ~ ldocs + lcpd | fe_id + year", s, "lcpd", f"size_q{q + 1}")
    out[f"size_q{q + 1}"]["jsize_range"] = [float(s.jsize.min()), float(s.jsize.max())]
pos["w"] = pos.docs_l1.clip(lower=1)
m = pf.feols("log_sjr ~ ldocs + lcpd | fe_id + year", data=pos, weights="w", vcov={"CRV1": "fe_id"}); t = m.tidy().loc["lcpd"]
out["weighted"] = {"est": float(t.Estimate), "se": float(t["Std. Error"]), "n": int(m._N)}; print("weighted", out["weighted"])
big = pos[pos.jsize >= 50]; est("log_sjr ~ ldocs + lcpd | fe_id + year", big, "lcpd", "size50"); out["size50"]["j"] = int(big.src_id.nunique())
est("log_sjr ~ ldocs + lcpd + lcpd_named + lcpd_ls | fe_id + year", pos, "lcpd_named", "gap_size")
for q in range(4):
    est("log_sjr ~ ldocs + lcpd + lcpd_named | fe_id + year", pos[pos.sq == q], "lcpd_named", f"gap_q{q + 1}")
json.dump(out, open(ROOT + "lc_size.json", "w"), indent=1)
