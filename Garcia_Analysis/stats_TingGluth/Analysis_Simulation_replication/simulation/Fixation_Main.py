# Combined fixation analyses (ES phase)
# - Duration (overall) ............................ FixDur_ES.csv + fixation_duration_ES_bars_SE.png
# - Probability (overall) ......................... Fix_Prop_ES.csv + fixation_location_probability_ES_bars_SE.png
# - Duration by Expected Value (E>S vs S>E) ....... fixation_duration_ES_by_value_condition_subjectmean.csv + fixation_duration_ES_by_value_BIG.png
# - Probability by Expected Value ................ FixLoc_Prob_ES_by_ExpectedValue.csv + fixation_location_probability_ES_by_value.png
# - Duration by EV × OV_2 ................... FixDur_ES_by_EV_OV.csv + per-OV plots
# - Probability by EV × OV_2 ................ FixLoc_Prob_ES_by_EV_OV.csv + per-OV plots
# All outputs go to ./fixtion_out

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

# ---------------------- config ----------------------
file = r"D:/Aberdeen_Uni_June24/cap/THESIS/Garcia_Analysis/data/data_sets/GarciaParticipants_Eye_Response_Feed_Allfix_addm_OV_Abs_CCT.csv"
excluded_subs = {1, 4, 5, 6, 14, 99}
OUT_DIR = "fixtion_out_OV"
os.makedirs(OUT_DIR, exist_ok=True)

FN_FIX_DUR_CSV  = "FixDur_ES.csv"
FN_FIX_DUR_PNG  = "fixation_duration_ES_bars_SE.png"
FN_FIX_PROP_CSV = "Fix_Prop_ES.csv"
FN_FIX_PROP_PNG = "fixation_location_probability_ES_bars_SE.png"
FN_DUR_VAL_CSV  = "fixation_duration_ES_by_value_condition_subjectmean.csv"
FN_DUR_VAL_PNG  = "fixation_duration_ES_by_value_BIG.png"
FN_PROP_VAL_CSV = "FixLoc_Prob_ES_by_ExpectedValue.csv"
FN_PROP_VAL_PNG = "fixation_location_probability_ES_by_value.png"
# EV × OV_2
FN_DUR_EV_OV_CSV  = "FixDur_ES_by_EV_OV.csv"
FN_PROP_EV_OV_CSV = "FixLoc_Prob_ES_by_EV_OV.csv"
# ----------------------------------------------------


# ---------------------- utils -----------------------
def to_numeric(x):
    try:
        if pd.isna(x):
            return np.nan
        return float(str(x).strip())
    except Exception:
        return np.nan

def loc_to_opt(loc_val):
    if pd.isna(loc_val):
        return np.nan
    v = int(round(loc_val))
    return "E" if v == 1 else ("S" if v == 2 else np.nan)

def sd1(x):
    return float(np.std(x, ddof=1)) if len(x) > 1 else np.nan

def add_log_columns(df, mean_col, sd_col, n_col):
  
    if "mean_log" in df.columns:
        out = df.copy()
        out["se_log"] = out["sd_log"] / np.sqrt(out[n_col].clip(lower=1))
        out["geom_mean_ms"] = np.exp(out["mu_log"])
        out["geom_sd_multiplier"] = np.exp(out["sigma_log"])
        return out
    else:
        # Fallback: if only raw mean/sd exist, do nothing extra
        return df

def compact_legend(ax, loc="upper right"):
    handles = [
        Patch(facecolor="blue", edgecolor="black", label="E (option)"),
        Patch(facecolor="red",  edgecolor="black", label="S (option)"),
        Patch(facecolor="white", edgecolor="black", label="Solid = E>S"),
        Patch(facecolor="white", edgecolor="black", hatch="//", label="Hatched = S>E"),
    ]
    ax.legend(handles=handles, loc=loc, frameon=False,
              prop={"size": 9}, labelspacing=0.3, borderpad=0.3,
              handletextpad=0.6, borderaxespad=0.3)
# ----------------------------------------------------


# ---------------------- load & preprocess ----------------------
df = pd.read_csv(file)

req = {"sub_id","phase",
       "FirstFixLoc","FirstFixDur","FinalFixLoc","FinalFixDur",
       "MiddleDominantLoc","MiddleDominantDur",
       "p1","p2","OV_2"}
missing = sorted(list(req - set(df.columns)))
if missing:
    raise ValueError(f"Missing required columns: {missing}")

df["sub_id"] = pd.to_numeric(df["sub_id"], errors="coerce").astype("Int64")
df["phase"]  = df["phase"].astype(str).str.upper()
for c in ["FirstFixLoc","FinalFixLoc","MiddleDominantLoc",
          "FirstFixDur","FinalFixDur","MiddleDominantDur",
          "p1","p2","OV_2"]:
    df[c] = df[c].map(to_numeric)

# ES-only, remove excluded subjects
df = df[(df["phase"] == "ES") & (~df["sub_id"].isin(excluded_subs))].copy()


def _find_each_middle_col(_df):
    # try common names (case-insensitive). Adjust if your actual name differs.
    candidates = [
        "eachMiddleFixDur", "eachMiddleFixation", "eachDwellMiddle",
        "EachMiddleFixDur", "each_middle_fix_dur", "each_middle_dur",
        "eachDwellMid", "eachDwellMiddleDur"
    ]
    clower = {c.lower(): c for c in _df.columns}
    for name in candidates:
        if name.lower() in clower:
            return clower[name.lower()]
    return None

def _fit_mu_sigma_from_series(ms_series):
    # filter valid positive durations, take natural logs, then mean/sd
    x = pd.to_numeric(ms_series, errors="coerce")
    x = x[(~x.isna()) & (x > 0)]
    if x.empty:
        return dict(n=0, mu=np.nan, sigma=np.nan, median_ms=np.nan)
    logs = np.log(x.values)
    mu = float(np.mean(logs))
    sigma = float(np.std(logs, ddof=1)) if logs.size > 1 else np.nan
    median_ms = float(np.exp(np.mean(logs)))  # geometric median (≈ typical ms)
    return dict(n=int(x.size), mu=mu, sigma=sigma, median_ms=median_ms)

first_stats = _fit_mu_sigma_from_series(df["FirstFixDur"])

EACH_MID_COL = _find_each_middle_col(df)
if EACH_MID_COL is None:
    print("Could not find an 'each middle' duration column. ")
    mid_stats = dict(n=0, mu=np.nan, sigma=np.nan, median_ms=np.nan)
else:
    mid_stats = _fit_mu_sigma_from_series(df[EACH_MID_COL])

# Save one tidy CSV for MATLAB
logparam_rows = [{
    "n_first": first_stats["n"],
    "mu_first": first_stats["mu"],
    "sigma_first": first_stats["sigma"],
    "median_first_ms": first_stats["median_ms"],
    "n_eachmiddle": mid_stats["n"],
    "mu_eachmiddle": mid_stats["mu"],
    "sigma_eachmiddle": mid_stats["sigma"],
    "median_eachmiddle_ms": mid_stats["median_ms"],
    "each_middle_colname": (EACH_MID_COL or "")
}]
logparam_csv = os.path.join(OUT_DIR, "FixDur_LogParams_ES_first_and_eachmiddle.csv")
pd.DataFrame(logparam_rows).to_csv(logparam_csv, index=False)
print("Saved:", os.path.abspath(logparam_csv))

#----------------------------------------------------------------------------------------

def _fit_mu_sigma(x_ms):
    x = pd.to_numeric(x_ms, errors="coerce")
    x = x[(~x.isna()) & (x > 0)]
    n = int(x.size)
    if n == 0:
        return dict(n=0, mu_log=np.nan, sigma_log=np.nan,
                    geom_median_ms=np.nan, geom_mean_ms=np.nan)
    logs = np.log(x.to_numpy())
    mu = float(np.mean(logs))
    sigma = float(np.std(logs, ddof=1)) if n > 1 else np.nan
    geom_median = float(np.exp(mu))
    geom_mean   = float(np.exp(mu + 0.5*(sigma**2))) if not np.isnan(sigma) else np.nan
    return dict(n=n, mu_log=mu, sigma_log=sigma,
                geom_median_ms=geom_median, geom_mean_ms=geom_mean)


rows_byov = []
for ov in [1, 2, 3]:
    sub = df[df["OV_2"] == ov]

    rows_byov.append({
        "OV_2": ov, "role": "First",
        **_fit_mu_sigma(sub["FirstFixDur"]),
        "source_col": "FirstFixDur"
    })

    if EACH_MID_COL is not None and EACH_MID_COL in sub.columns:
        rows_byov.append({
            "OV_2": ov, "role": "EachMiddle",
            **_fit_mu_sigma(sub[EACH_MID_COL]),
            "source_col": EACH_MID_COL
        })
    else:
        rows_byov.append({
            "OV_2": ov, "role": "EachMiddle",
            "n": 0, "mu_log": np.nan, "sigma_log": np.nan,
            "geom_median_ms": np.nan, "geom_mean_ms": np.nan,
            "source_col": EACH_MID_COL or ""
        })

FN_LOG_BYOV = os.path.join(OUT_DIR, "FixDur_LogParams_ES_first_and_eachmiddle_by_OV.csv")
pd.DataFrame(rows_byov).to_csv(FN_LOG_BYOV, index=False)
print("Saved:", os.path.abspath(FN_LOG_BYOV))

# Ensure p1 and p2 are numeric and drop problematic rows
df["p1"] = pd.to_numeric(df["p1"], errors="coerce")
df["p2"] = pd.to_numeric(df["p2"], errors="coerce")
df = df.dropna(subset=["p1", "p2"]).copy()

TIE_TOL = 1e-3
diff = df["p1"] - df["p2"]
cond3 = np.where(diff >  TIE_TOL, "E>S",
         np.where(diff < -TIE_TOL, "S>E", "E=S"))
df["value_cond3"] = pd.Categorical(cond3, ["E>S","E=S","S>E"], ordered=True)

cond = np.select(
    [df["p1"] > df["p2"], df["p2"] > df["p1"]],
    ["E>S", "S>E"],
    default=None
)

df["value_cond"] = pd.Categorical(cond, ["E>S","S>E"], ordered=True)

EACH_MID_COL = _find_each_middle_col(df)
if EACH_MID_COL is None:
    raise ValueError(
        "Each-middle fixation duration column not found. "
        "Add it or set `EACH_MID_COL = 'your_column_name'` before make_dur_long()."
    )

def make_dur_long(d):
    first  = pd.DataFrame({"sub_id": d["sub_id"], "role":"First",
                           "opt_type": d["FirstFixLoc"].map(loc_to_opt),
                           "dur_ms": d["FirstFixDur"],
                           "value_cond": d["value_cond"],
                           "value_cond3": d["value_cond3"],
                           "OV_2": d["OV_2"]})
    middle = pd.DataFrame({"sub_id": d["sub_id"], "role":"Middle",
                           "opt_type": d["MiddleDominantLoc"].map(loc_to_opt),
                           "dur_ms": d[EACH_MID_COL],
                           "value_cond": d["value_cond"],
                           "value_cond3": d["value_cond3"],
                           "OV_2": d["OV_2"]})
    final  = pd.DataFrame({"sub_id": d["sub_id"], "role":"Final",
                           "opt_type": d["FinalFixLoc"].map(loc_to_opt),
                           "dur_ms": d["FinalFixDur"],
                           "value_cond": d["value_cond"],
                           "value_cond3": d["value_cond3"],
                           "OV_2": d["OV_2"]})

    out = pd.concat([first, middle, final], ignore_index=True)
    out = out[out["opt_type"].isin(["E","S"]) & (~out["dur_ms"].isna()) & (out["dur_ms"] > 0)].copy()
    out["role"] = pd.Categorical(out["role"], ["First","Middle","Final"], ordered=True)
    out["opt_type"] = pd.Categorical(out["opt_type"], ["E","S"], ordered=True)
    out["OV_2"] = out["OV_2"].astype("Int64")
    out["log_dur"] = np.log(out["dur_ms"])
    return out


def make_loc_long(d):
    first  = pd.DataFrame({"sub_id": d["sub_id"], "role":"First",
                           "opt_type": d["FirstFixLoc"].map(loc_to_opt),
                           "value_cond": d["value_cond"],
                           "value_cond3": d["value_cond3"],  
                           "OV_2": d["OV_2"]})
    middle = pd.DataFrame({"sub_id": d["sub_id"], "role":"Middle",
                           "opt_type": d["MiddleDominantLoc"].map(loc_to_opt),
                           "value_cond": d["value_cond"],
                           "value_cond3": d["value_cond3"],   
                           "OV_2": d["OV_2"]})
    final  = pd.DataFrame({"sub_id": d["sub_id"], "role":"Final",
                           "opt_type": d["FinalFixLoc"].map(loc_to_opt),
                           "value_cond": d["value_cond"],
                           "value_cond3": d["value_cond3"],  
                           "OV_2": d["OV_2"]})
    out = pd.concat([first, middle, final], ignore_index=True)
    out = out[out["opt_type"].isin(["E","S"])].copy()
    out["role"] = pd.Categorical(out["role"], ["First","Middle","Final"], ordered=True)
    out["opt_type"] = pd.Categorical(out["opt_type"], ["E","S"], ordered=True)
    out["OV_2"] = out["OV_2"].astype("Int64")
    return out

dur_long = make_dur_long(df)
loc_long = make_loc_long(df)

roles = ["First","Middle","Final"]
opt_types = ["E","S"]
role_idx = {r:i for i,r in enumerate(roles)}
# ---------------------------------------------------------------


# Fixation Durations (overall) ----------------------
# per-subject (raw & log)
subj_dur = (dur_long
    .groupby(["sub_id","role","opt_type"], observed=True, as_index=False)
    .agg(mean_ms=("dur_ms","mean"),
         median_ms=("dur_ms","median"),
         mean_log=("log_dur","mean"),
         sd_log=("log_dur","std"),
         n_trials=("dur_ms","size"))
)

# across-subject (raw) + (log summaries)
dur_all = (subj_dur
    .groupby(["role","opt_type"], observed=True, as_index=False)
    .agg(n_subjects=("sub_id","nunique"),
         n_trials_total=("n_trials","sum"),
         mean_ms=("mean_ms","mean"),
         sd_ms=("mean_ms", sd1),
         median_ms=("median_ms","median"),
         mu_log=("mean_log","mean"),
         sigma_log=("mean_log", sd1))   # SD across subjects of per-subject mean(log)
)
dur_all["se_ms"] = dur_all["sd_ms"] / np.sqrt(dur_all["n_subjects"].clip(lower=1))
dur_all["se_log"] = dur_all["sigma_log"] / np.sqrt(dur_all["n_subjects"].clip(lower=1))
dur_all["geom_mean_ms"] = np.exp(dur_all["mu_log"])
dur_all["geom_sd_multiplier"] = np.exp(dur_all["sigma_log"])

dur_all = dur_all.sort_values(["role","opt_type"])
dur_all.to_csv(os.path.join(OUT_DIR, FN_FIX_DUR_CSV), index=False)
print("Saved:", os.path.abspath(os.path.join(OUT_DIR, FN_FIX_DUR_CSV)))

# plot
mean_grid = np.full((len(roles), 2), np.nan); se_grid = np.full_like(mean_grid, np.nan, dtype=float)
for _, row in dur_all.iterrows():
    i = role_idx[row["role"]]; j = opt_types.index(row["opt_type"])
    mean_grid[i, j] = row["mean_ms"]; se_grid[i, j] = row["se_ms"]

x = np.arange(len(roles)); width = 0.35
fig, ax = plt.subplots(figsize=(8,5), dpi=150)
ax.bar(x - width/2, mean_grid[:,0], width, label="E", color="blue", yerr=np.nan_to_num(se_grid[:,0]), capsize=4)
ax.bar(x + width/2, mean_grid[:,1], width, label="S", color="red",  yerr=np.nan_to_num(se_grid[:,1]), capsize=4)
ax.set_ylabel("Fixation Duration (ms)")
ax.set_xticks(x, roles)
ax.set_title("Fixation Durations (mean ± SE)")
ax.legend(loc="upper left", ncols=2)
ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
fig.tight_layout()
fig.savefig(os.path.join(OUT_DIR, FN_FIX_DUR_PNG))
plt.close(fig)
print("Saved:", os.path.abspath(os.path.join(OUT_DIR, FN_FIX_DUR_PNG)))


# Fixation Location Probability (overall) ----------------------
counts = (loc_long
          .groupby(["sub_id","role","opt_type"])
          .size()
          .unstack(fill_value=0)
          .reindex(columns=["E","S"], fill_value=0)
          .reset_index())
counts["total"] = counts["E"] + counts["S"]
counts = counts[counts["total"] > 0].copy()
counts["prob_E"] = counts["E"] / counts["total"]
counts["prob_S"] = counts["S"] / counts["total"]

subj_probs = counts.melt(
    id_vars=["sub_id","role"], value_vars=["prob_E","prob_S"],
    var_name="opt_type", value_name="prob"
)
subj_probs["opt_type"] = subj_probs["opt_type"].map({"prob_E":"E","prob_S":"S"})
subj_probs["opt_type"] = pd.Categorical(subj_probs["opt_type"], ["E","S"], ordered=True)

prop_all = (subj_probs
    .groupby(["role","opt_type"], observed=True, as_index=False)
    .agg(n_subjects=("sub_id","nunique"),
         mean_prob=("prob","mean"),
         sd_prob=("prob", sd1),
         median_prob=("prob","median"))
)
prop_all["se_prob"]      = prop_all["sd_prob"] / np.sqrt(prop_all["n_subjects"].clip(lower=1))
prop_all["mean_percent"] = prop_all["mean_prob"] * 100.0
prop_all["se_percent"]   = prop_all["se_prob"]   * 100.0
prop_all["median_percent"] = prop_all["median_prob"] * 100.0

prop_all = prop_all.sort_values(["role","opt_type"])
prop_all.to_csv(os.path.join(OUT_DIR, FN_FIX_PROP_CSV), index=False)
print("Saved:", os.path.abspath(os.path.join(OUT_DIR, FN_FIX_PROP_CSV)))

# plot
mean_grid = np.full((len(roles), 2), np.nan); se_grid = np.full_like(mean_grid, np.nan, dtype=float)
for _, row in prop_all.iterrows():
    i = role_idx[row["role"]]; j = opt_types.index(row["opt_type"])
    mean_grid[i, j] = row["mean_percent"]; se_grid[i, j] = row["se_percent"]

fig, ax = plt.subplots(figsize=(8,5), dpi=150)
ax.bar(x - width/2, mean_grid[:,0], width, label="E", color="blue", yerr=np.nan_to_num(se_grid[:,0]), capsize=4)
ax.bar(x + width/2, mean_grid[:,1], width, label="S", color="red",  yerr=np.nan_to_num(se_grid[:,1]), capsize=4)
ax.set_ylabel("Probability (%)")
ax.set_xticks(x, roles); ax.set_ylim(0, 100)
ax.set_title("Probability Fixation Falls on E vs S (mean ± SE)")
ax.legend(loc="upper center", ncols=2)
ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
fig.tight_layout()
fig.savefig(os.path.join(OUT_DIR, FN_FIX_PROP_PNG))
plt.close(fig)
print("Saved:", os.path.abspath(os.path.join(OUT_DIR, FN_FIX_PROP_PNG)))


# Fixation Durations by Expected Value ----------------------
dur_val = dur_long.dropna(subset=["value_cond"]).copy()

# per-subject (raw+log), grouped by EV condition
subj_val = (dur_val
    .groupby(["sub_id","role","value_cond","opt_type"], observed=True, as_index=False)
    .agg(mean_ms=("dur_ms","mean"),
         median_ms=("dur_ms","median"),
         mean_log=("log_dur","mean"),
         sd_log=("log_dur","std"),
         n_trials=("dur_ms","size"))
)

# across-subject
dur_ev = (subj_val
    .groupby(["role","value_cond","opt_type"], observed=True, as_index=False)
    .agg(n_subjects=("sub_id","nunique"),
         n_trials_total=("n_trials","sum"),
         mean_ms=("mean_ms","mean"),
         sd_ms=("mean_ms", sd1),
         median_ms=("median_ms","median"),
         mu_log=("mean_log","mean"),
         sigma_log=("mean_log", sd1))
    .sort_values(["role","value_cond","opt_type"])
)
dur_ev["se_ms"]  = dur_ev["sd_ms"]  / np.sqrt(dur_ev["n_subjects"].clip(lower=1))
dur_ev["se_log"] = dur_ev["sigma_log"] / np.sqrt(dur_ev["n_subjects"].clip(lower=1))
dur_ev["geom_mean_ms"] = np.exp(dur_ev["mu_log"])
dur_ev["geom_sd_multiplier"] = np.exp(dur_ev["sigma_log"])

dur_ev.to_csv(os.path.join(OUT_DIR, FN_DUR_VAL_CSV), index=False)
print("Saved:", os.path.abspath(os.path.join(OUT_DIR, FN_DUR_VAL_CSV)))

# big plot 
mapping = [("E>S","E"), ("E>S","S"), ("S>E","E"), ("S>E","S")]
mean_mat = np.full((len(roles), 4), np.nan); se_mat = np.full_like(mean_mat, np.nan, dtype=float)
for i, r in enumerate(roles):
    for k, (c,o) in enumerate(mapping):
        row = dur_ev[(dur_ev["role"]==r) & (dur_ev["value_cond"]==c) & (dur_ev["opt_type"]==o)]
        if len(row) == 1:
            mean_mat[i, k] = row["mean_ms"].values[0]
            se_mat[i, k]   = row["se_ms"].values[0]

x = np.arange(len(roles)); width_small = 0.18
pos = [x - 1.5*width_small, x - 0.5*width_small, x + 0.5*width_small, x + 1.5*width_small]

fig, ax = plt.subplots(figsize=(10,5.8), dpi=150)
ax.bar(pos[0], mean_mat[:,0], width_small, yerr=np.nan_to_num(se_mat[:,0]), capsize=4, color="blue")
ax.bar(pos[1], mean_mat[:,1], width_small, yerr=np.nan_to_num(se_mat[:,1]), capsize=4, color="red")
ax.bar(pos[2], mean_mat[:,2], width_small, yerr=np.nan_to_num(se_mat[:,2]), capsize=4, color="blue", hatch="//")
ax.bar(pos[3], mean_mat[:,3], width_small, yerr=np.nan_to_num(se_mat[:,3]), capsize=4, color="red",  hatch="//")

ax.set_xticks(x, roles)
ax.set_ylabel("Fixation Duration (ms)")
ax.set_title("Fixation Duration by Value Dominance and Option (mean ± SE)", pad=10)
ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
compact_legend(ax, loc="upper left")

fig.tight_layout()
fig.savefig(os.path.join(OUT_DIR, FN_DUR_VAL_PNG))
plt.close(fig)
print("Saved:", os.path.abspath(os.path.join(OUT_DIR, FN_DUR_VAL_PNG)))


# Fixation Location Probability by Expected Value ----------------------
loc_val = loc_long.dropna(subset=["value_cond3"]).copy()

counts_val = (loc_val.groupby(["sub_id","role","value_cond3","opt_type"])
              .size().unstack(fill_value=0).reindex(columns=["E","S"], fill_value=0).reset_index())

counts_val["total"]  = counts_val["E"] + counts_val["S"]
counts_val = counts_val[counts_val["total"] > 0].copy()
counts_val["prob_E"] = counts_val["E"] / counts_val["total"]
counts_val["prob_S"] = counts_val["S"] / counts_val["total"]

subj_probs_val = counts_val.melt(
    id_vars=["sub_id","role","value_cond3"], value_vars=["prob_E","prob_S"],
    var_name="opt_type", value_name="prob"
)
subj_probs_val["opt_type"] = subj_probs_val["opt_type"].map({"prob_E":"E","prob_S":"S"})
subj_probs_val["opt_type"] = pd.Categorical(subj_probs_val["opt_type"], ["E","S"], ordered=True)

prop_ev = (subj_probs_val
    .groupby(["role","value_cond3","opt_type"], observed=True, as_index=False)
    .agg(n_subjects=("sub_id","nunique"),
         mean_prob=("prob","mean"),
         sd_prob=("prob", sd1),
         median_prob=("prob","median"))
    .sort_values(["role","value_cond3","opt_type"])
)

prop_ev["se_prob"]      = prop_ev["sd_prob"] / np.sqrt(prop_ev["n_subjects"].clip(lower=1))
prop_ev["mean_percent"] = prop_ev["mean_prob"] * 100.0
prop_ev["se_percent"]   = prop_ev["se_prob"]   * 100.0
prop_ev["median_percent"] = prop_ev["median_prob"] * 100.0

prop_ev.to_csv(os.path.join(OUT_DIR, FN_PROP_VAL_CSV), index=False)
print("Saved:", os.path.abspath(os.path.join(OUT_DIR, FN_PROP_VAL_CSV)))

# big plot
mean_mat = np.full((len(roles), 4), np.nan); se_mat = np.full_like(mean_mat, np.nan, dtype=float)
for i, r in enumerate(roles):
    for k, (c,o) in enumerate(mapping):   
        row = prop_ev[(prop_ev["role"]==r) & (prop_ev["value_cond3"]==c) & (prop_ev["opt_type"]==o)]
        if len(row) == 1:
            mean_mat[i, k] = row["mean_percent"].values[0]
            se_mat[i, k]   = row["se_percent"].values[0]

fig, ax = plt.subplots(figsize=(10,5.8), dpi=150)
ax.bar(pos[0], mean_mat[:,0], width_small, yerr=np.nan_to_num(se_mat[:,0]), capsize=4, color="blue")
ax.bar(pos[1], mean_mat[:,1], width_small, yerr=np.nan_to_num(se_mat[:,1]), capsize=4, color="red")
ax.bar(pos[2], mean_mat[:,2], width_small, yerr=np.nan_to_num(se_mat[:,2]), capsize=4, color="blue", hatch="//")
ax.bar(pos[3], mean_mat[:,3], width_small, yerr=np.nan_to_num(se_mat[:,3]), capsize=4, color="red",  hatch="//")

ax.set_xticks(x, roles)
ax.set_ylabel("Probability (%)")
ax.set_ylim(0, 100)
ax.set_title("P(Fixation on E vs S) by Expected Value and Option (mean ± SE)", pad=10)
ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
compact_legend(ax, loc="upper right")

fig.tight_layout()
fig.savefig(os.path.join(OUT_DIR, FN_PROP_VAL_PNG))
plt.close(fig)
print("Saved:", os.path.abspath(os.path.join(OUT_DIR, FN_PROP_VAL_PNG)))


# Duration by EV nd OV_2 (with log summaries) 
dur_ev_ov = dur_long.dropna(subset=["value_cond","OV_2"]).copy()

subj_ev_ov = (dur_ev_ov
    .groupby(["sub_id","role","value_cond","OV_2","opt_type"], observed=True, as_index=False)
    .agg(mean_ms=("dur_ms","mean"),
         median_ms=("dur_ms","median"),
         mean_log=("log_dur","mean"),
         sd_log=("log_dur","std"),
         n_trials=("dur_ms","size"))
)

dur_ev_ov_sum = (subj_ev_ov
    .groupby(["role","value_cond","OV_2","opt_type"], observed=True, as_index=False)
    .agg(n_subjects=("sub_id","nunique"),
         n_trials_total=("n_trials","sum"),
         mean_ms=("mean_ms","mean"),
         sd_ms=("mean_ms", sd1),
         median_ms=("median_ms","median"),
         mu_log=("mean_log","mean"),
         sigma_log=("mean_log", sd1))
    .sort_values(["role","value_cond","OV_2","opt_type"])
)
dur_ev_ov_sum["se_ms"]  = dur_ev_ov_sum["sd_ms"]  / np.sqrt(dur_ev_ov_sum["n_subjects"].clip(lower=1))
dur_ev_ov_sum["se_log"] = dur_ev_ov_sum["sigma_log"] / np.sqrt(dur_ev_ov_sum["n_subjects"].clip(lower=1))
dur_ev_ov_sum["geom_mean_ms"] = np.exp(dur_ev_ov_sum["mu_log"])
dur_ev_ov_sum["geom_sd_multiplier"] = np.exp(dur_ev_ov_sum["sigma_log"])

dur_ev_ov_sum.to_csv(os.path.join(OUT_DIR, FN_DUR_EV_OV_CSV), index=False)
print("Saved:", os.path.abspath(os.path.join(OUT_DIR, FN_DUR_EV_OV_CSV)))

# per-OV plots (one figure per OV_2 level)
for ov in [1,2,3]:
    tmp = dur_ev_ov_sum[dur_ev_ov_sum["OV_2"] == ov]
    if tmp.empty: 
        continue
    mean_mat = np.full((len(roles), 4), np.nan); se_mat = np.full_like(mean_mat, np.nan, dtype=float)
    for i, r in enumerate(roles):
        for k, (c,o) in enumerate(mapping):
            row = tmp[(tmp["role"]==r) & (tmp["value_cond"]==c) & (tmp["opt_type"]==o)]
            if len(row) == 1:
                mean_mat[i, k] = row["mean_ms"].values[0]
                se_mat[i, k]   = row["se_ms"].values[0]
    fig, ax = plt.subplots(figsize=(10,5.8), dpi=150)
    ax.bar(pos[0], mean_mat[:,0], width_small, yerr=np.nan_to_num(se_mat[:,0]), capsize=4, color="blue")
    ax.bar(pos[1], mean_mat[:,1], width_small, yerr=np.nan_to_num(se_mat[:,1]), capsize=4, color="red")
    ax.bar(pos[2], mean_mat[:,2], width_small, yerr=np.nan_to_num(se_mat[:,2]), capsize=4, color="blue", hatch="//")
    ax.bar(pos[3], mean_mat[:,3], width_small, yerr=np.nan_to_num(se_mat[:,3]), capsize=4, color="red",  hatch="//")
    ax.set_xticks(x, roles)
    ax.set_ylabel("Fixation Duration (ms)")
    ax.set_title(f"Fixation Duration by EV and Option (OV={ov}; mean ± SE)", pad=10)
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    compact_legend(ax, loc="upper left")
    fig.tight_layout()
    fig.savefig(os.path.join(OUT_DIR, f"fixation_duration_ES_by_value_by_OV{ov}.png"))
    plt.close(fig)
    print("Saved:", os.path.abspath(os.path.join(OUT_DIR, f"fixation_duration_ES_by_value_by_OV{ov}.png")))


# Probability by EV × OV_2 ----------------------
loc_ev_ov = loc_long.dropna(subset=["value_cond3","OV_2"]).copy()

counts_ev_ov = (loc_ev_ov.groupby(["sub_id","role","value_cond3","OV_2","opt_type"])
                .size().unstack(fill_value=0).reindex(columns=["E","S"], fill_value=0).reset_index())
counts_ev_ov["total"]  = counts_ev_ov["E"] + counts_ev_ov["S"]
counts_ev_ov = counts_ev_ov[counts_ev_ov["total"] > 0].copy()
counts_ev_ov["prob_E"] = counts_ev_ov["E"] / counts_ev_ov["total"]
counts_ev_ov["prob_S"] = counts_ev_ov["S"] / counts_ev_ov["total"]

subj_probs_ev_ov = counts_ev_ov.melt(
    id_vars=["sub_id","role","value_cond3","OV_2"], value_vars=["prob_E","prob_S"],
    var_name="opt_type", value_name="prob"
)
subj_probs_ev_ov["opt_type"] = subj_probs_ev_ov["opt_type"].map({"prob_E":"E","prob_S":"S"})
subj_probs_ev_ov["opt_type"] = pd.Categorical(subj_probs_ev_ov["opt_type"], ["E","S"], ordered=True)

prop_ev_ov = (subj_probs_ev_ov
    .groupby(["role","value_cond3","OV_2","opt_type"], observed=True, as_index=False)
    .agg(n_subjects=("sub_id","nunique"),
         mean_prob=("prob","mean"),
         sd_prob=("prob", sd1),
         median_prob=("prob","median"))
    .sort_values(["role","value_cond3","OV_2","opt_type"])   # <-- here
)

prop_ev_ov["se_prob"]      = prop_ev_ov["sd_prob"] / np.sqrt(prop_ev_ov["n_subjects"].clip(lower=1))
prop_ev_ov["mean_percent"] = prop_ev_ov["mean_prob"] * 100.0
prop_ev_ov["se_percent"]   = prop_ev_ov["se_prob"]   * 100.0
prop_ev_ov["median_percent"] = prop_ev_ov["median_prob"] * 100.0

prop_ev_ov.to_csv(os.path.join(OUT_DIR, FN_PROP_EV_OV_CSV), index=False)
print("Saved:", os.path.abspath(os.path.join(OUT_DIR, FN_PROP_EV_OV_CSV)))

# per-OV plots
for ov in [1,2,3]:
    tmp = prop_ev_ov[prop_ev_ov["OV_2"] == ov]
    # (optional) explicitly drop ties from the plotting DataFrame:
    tmp = tmp[tmp["value_cond3"].isin(["E>S","S>E"])]
    if tmp.empty:
        continue

    mean_mat = np.full((len(roles), 4), np.nan)
    se_mat   = np.full_like(mean_mat, np.nan, dtype=float)

    for i, r in enumerate(roles):
        for k, (c,o) in enumerate(mapping):
            row = tmp[(tmp["role"]==r) & (tmp["value_cond3"]==c) & (tmp["opt_type"]==o)]
            if len(row) == 1:
                mean_mat[i, k] = row["mean_percent"].values[0]
                se_mat[i, k]   = row["se_percent"].values[0]

    fig, ax = plt.subplots(figsize=(10,5.8), dpi=150)
    ax.bar(pos[0], mean_mat[:,0], width_small, yerr=np.nan_to_num(se_mat[:,0]), capsize=4, color="blue")
    ax.bar(pos[1], mean_mat[:,1], width_small, yerr=np.nan_to_num(se_mat[:,1]), capsize=4, color="red")
    ax.bar(pos[2], mean_mat[:,2], width_small, yerr=np.nan_to_num(se_mat[:,2]), capsize=4, color="blue", hatch="//")
    ax.bar(pos[3], mean_mat[:,3], width_small, yerr=np.nan_to_num(se_mat[:,3]), capsize=4, color="red",  hatch="//")
    ax.set_xticks(x, roles)
    ax.set_ylabel("Probability (%)")
    ax.set_ylim(0, 100)
    ax.set_title(f"P(Fixation on E vs S) by EV and Option (OV={ov}; mean ± SE)", pad=10)
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    compact_legend(ax, loc="upper right")
    fig.tight_layout()
    fig.savefig(os.path.join(OUT_DIR, f"fixation_location_probability_ES_by_value_by_OV{ov}.png"))
    plt.close(fig)
    print("Saved:", os.path.abspath(os.path.join(OUT_DIR, f"fixation_location_probability_ES_by_value_by_OV{ov}.png")))

print("\nDone. Files in:", os.path.abspath(OUT_DIR))


####----------------------------------------------------------------------------------------------------------------------

# Duration and Probability by EV and OV_VD_label

print("\nComputing Duration and Probability by EV and OV_VD_label")

# Make sure OV_VD_label exists (create if not)
if "OV_VD_label" not in df.columns:
    df = df.dropna(subset=["OV_2", "VD_2"]).copy()
    df["OV_2"] = df["OV_2"].astype(int)
    df["VD_2"] = df["VD_2"].astype(int)
    df["OV_VD_label"] = (df["OV_2"].astype(str) + df["VD_2"].astype(str)).astype(int)

# --- DURATIONS ---
dur_ev_ovvd = dur_long.dropna(subset=["value_cond", "OV_VD_label"]).copy()

subj_ev_ovvd = (
    dur_ev_ovvd.groupby(["sub_id", "role", "value_cond", "OV_VD_label", "opt_type"],
                        observed=True, as_index=False)
    .agg(mean_ms=("dur_ms", "mean"),
         median_ms=("dur_ms", "median"),
         mean_log=("log_dur", "mean"),
         sd_log=("log_dur", "std"),
         n_trials=("dur_ms", "size"))
)

dur_ev_ovvd_sum = (
    subj_ev_ovvd.groupby(["role", "value_cond", "OV_VD_label", "opt_type"],
                         observed=True, as_index=False)
    .agg(n_subjects=("sub_id", "nunique"),
         n_trials_total=("n_trials", "sum"),
         mean_ms=("mean_ms", "mean"),
         sd_ms=("mean_ms", sd1),
         median_ms=("median_ms", "median"),
         mu_log=("mean_log", "mean"),
         sigma_log=("mean_log", sd1))
    .sort_values(["role", "value_cond", "OV_VD_label", "opt_type"])
)

dur_ev_ovvd_sum["se_ms"]  = dur_ev_ovvd_sum["sd_ms"]  / np.sqrt(dur_ev_ovvd_sum["n_subjects"].clip(lower=1))
dur_ev_ovvd_sum["se_log"] = dur_ev_ovvd_sum["sigma_log"] / np.sqrt(dur_ev_ovvd_sum["n_subjects"].clip(lower=1))
dur_ev_ovvd_sum["geom_mean_ms"] = np.exp(dur_ev_ovvd_sum["mu_log"])
dur_ev_ovvd_sum["geom_sd_multiplier"] = np.exp(dur_ev_ovvd_sum["sigma_log"])

FN_DUR_EV_OVVD_CSV = os.path.join(OUT_DIR, "FixDur_ES_by_EV_OVVDlabel.csv")
dur_ev_ovvd_sum.to_csv(FN_DUR_EV_OVVD_CSV, index=False)
print("Saved:", os.path.abspath(FN_DUR_EV_OVVD_CSV))


# PROBABILITIES 
loc_ev_ovvd = loc_long.dropna(subset=["value_cond3", "OV_VD_label"]).copy()

counts_ev_ovvd = (
    loc_ev_ovvd.groupby(["sub_id", "role", "value_cond3", "OV_VD_label", "opt_type"])
    .size().unstack(fill_value=0).reindex(columns=["E", "S"], fill_value=0).reset_index()
)
counts_ev_ovvd["total"] = counts_ev_ovvd["E"] + counts_ev_ovvd["S"]
counts_ev_ovvd = counts_ev_ovvd[counts_ev_ovvd["total"] > 0].copy()
counts_ev_ovvd["prob_E"] = counts_ev_ovvd["E"] / counts_ev_ovvd["total"]
counts_ev_ovvd["prob_S"] = counts_ev_ovvd["S"] / counts_ev_ovvd["total"]

subj_probs_ev_ovvd = counts_ev_ovvd.melt(
    id_vars=["sub_id", "role", "value_cond3", "OV_VD_label"],
    value_vars=["prob_E", "prob_S"],
    var_name="opt_type", value_name="prob"
)
subj_probs_ev_ovvd["opt_type"] = subj_probs_ev_ovvd["opt_type"].map({"prob_E": "E", "prob_S": "S"})
subj_probs_ev_ovvd["opt_type"] = pd.Categorical(subj_probs_ev_ovvd["opt_type"], ["E", "S"], ordered=True)

prop_ev_ovvd = (
    subj_probs_ev_ovvd.groupby(["role", "value_cond3", "OV_VD_label", "opt_type"],
                               observed=True, as_index=False)
    .agg(n_subjects=("sub_id", "nunique"),
         mean_prob=("prob", "mean"),
         sd_prob=("prob", sd1),
         median_prob=("prob", "median"))
    .sort_values(["role", "value_cond3", "OV_VD_label", "opt_type"])
)
prop_ev_ovvd["se_prob"]      = prop_ev_ovvd["sd_prob"] / np.sqrt(prop_ev_ovvd["n_subjects"].clip(lower=1))
prop_ev_ovvd["mean_percent"] = prop_ev_ovvd["mean_prob"] * 100.0
prop_ev_ovvd["se_percent"]   = prop_ev_ovvd["se_prob"]   * 100.0
prop_ev_ovvd["median_percent"] = prop_ev_ovvd["median_prob"] * 100.0

FN_PROP_EV_OVVD_CSV = os.path.join(OUT_DIR, "FixLoc_Prob_ES_by_EV_OVVDlabel.csv")
prop_ev_ovvd.to_csv(FN_PROP_EV_OVVD_CSV, index=False)
print("Saved:", os.path.abspath(FN_PROP_EV_OVVD_CSV))

print("\n Saved new OV_VD_label-based CSVs.")
