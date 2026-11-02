"""4_fe_checks.py
Two checks estimated with pyfixest. Run after 1_build_panel.R and 2_models.R.
  (a) primary field x year fixed effects (Table 9, Section 3.6)
  (b) density moderation (M8, strict) including journals indexed from 1996 or earlier
  (c) composition checks: year slopes on journals observed in every year; field deviations with early journals
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

# ---------- (b) density moderation including early journals
def prep(path):
    d = pd.read_csv(path); fv = [c for c in d.columns if c.startswith("field_")]
    d["fe_id"] = d.src_id.astype(str); d["ldocs"] = np.log1p(d.docs_l1)
    d = d[d.cpd_l1 > 0].copy(); d["lcpd"] = np.log(d.cpd_l1)
    G = d[fv].values.astype(float); nf = np.maximum(G.sum(1), 1)
    dens = np.zeros(len(d)); refd = np.zeros(len(d)); d["lrpd"] = np.log(d.rpd_l1.clip(lower=0.5))
    for j, f in enumerate(fv):
        m = d[f] == 1
        dens += G[:, j] * d.year.map(d[m].groupby("year").lcpd.mean()).fillna(0).values
        refd += G[:, j] * d.year.map(d[m].groupby("year").lrpd.mean()).fillna(0).values
    d["dens_c"] = dens / nf; d["dens_c"] -= d.dens_c.mean(); d["ref_c"] = refd / nf; d["ref_c"] -= d.ref_c.mean()
    d["lsq"] = (d.lcpd - d.lcpd.mean())**2
    # cluster id for field-clustered SEs: first of the journal's fields in field_vars order
    d["pfield"] = np.argmax(G, axis=1).astype(str); d["fy"] = d.pfield + "_" + d.year.astype(str)
    yrs = sorted(d.year.unique())[1:]
    for y in yrs: d[f"ly{y}"] = d.lcpd * (d.year == y)
    return d, fv, yrs
full, fvf, yrsf = prep(ROOT + "full_sample.csv")
fif = " + ".join(f"lcpd:{f}" for f in fvf); lyf = " + ".join(f"ly{y}" for y in yrsf)
out = {}
for dv in ["dens_c", "ref_c"]:
    m = pf.feols(f"log_sjr ~ ldocs + lcpd + {dv} + lcpd:{dv} + lsq + {fif} + {lyf} | fe_id + year", data=full, vcov={"CRV1": "fe_id"})
    t = m.tidy().loc[f"lcpd:{dv}"]; out[dv] = {"est": t.Estimate, "se_journal": t["Std. Error"]}
    m.vcov({"CRV1": "pfield"}); t = m.tidy().loc[f"lcpd:{dv}"]; out[dv]["se_field"] = t["Std. Error"]; out[dv]["p_field"] = t["Pr(>|t|)"]
    print(dv, out[dv], flush=True); del m; gc.collect()
json.dump(out, open(ROOT + "lc_full_density.json", "w"), indent=1)

# ---------- (c) composition checks (Section 3.6)
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
