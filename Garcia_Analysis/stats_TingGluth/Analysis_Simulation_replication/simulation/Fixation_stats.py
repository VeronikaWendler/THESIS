# Get Fixation Stats

#libs
import os
import numpy as np
import pandas as pd

file = r"D:/Aberdeen_Uni_June24/cap/THESIS/Garcia_Analysis/data/data_sets/GarciaParticipants_Eye_Response_Feed_Allfix_addm_OV_Abs_CCT.csv"
excluded_subs = {1, 4, 5, 6, 14, 99}   # for EE: 6, 14, 20, 26, 2, 9, 18

PHASE = "ES"  # determine phase
OUT_DIR = "fixation_out_ES"
os.makedirs(OUT_DIR, exist_ok=True)

FN_DUR_EV_OV  = f"FixDur_{PHASE}_by_EV_OV.csv"
FN_PROP_EV_OV = f"FixLoc_Prob_{PHASE}_by_EV_OV.csv"

# helpers ----------------------------------------------------------------------
def to_numeric(x):
    try:
        if pd.isna(x):
            return np.nan
        return float(str(x).strip())
    except Exception:
        return np.nan

def loc_to_side(loc_val):
    if pd.isna(loc_val):
        return np.nan
    v = int(round(loc_val))
    return "left" if v == 1 else ("right" if v == 2 else np.nan)

def sd1(x):
    x = np.asarray(x, dtype=float)
    return float(np.std(x, ddof=1)) if x.size > 1 else np.nan

# ---------------- prep -------------------------------------------------------

df = pd.read_csv(file)

required = [
    "sub_id","phase","OV_2","p1","p2",
    "FirstFixLoc","FirstFixDur",
    "MiddleDominantLoc","eachMiddleFixDur",
    "FinalFixLoc","FinalFixDur",
]

missing = [c for c in required if c not in df.columns]
if missing:
    raise ValueError(f"Missing required columns: {missing}")

df["sub_id"] = pd.to_numeric(df["sub_id"], errors="coerce").astype("Int64")
df["phase"]  = df["phase"].astype(str).str.upper()

for c in [
    "OV_2","p1","p2",
    "FirstFixLoc","FirstFixDur",
    "MiddleDominantLoc","eachMiddleFixDur",
    "FinalFixLoc","FinalFixDur"
]:
    df[c] = df[c].map(to_numeric)

df = df[(df["phase"] == PHASE) & (~df["sub_id"].isin(excluded_subs))].copy()

df = df.dropna(subset=["OV_2","p1","p2"]).copy()
df["OV_2"] = df["OV_2"].astype("Int64")

# ---------------- expected value check 
diff = df["p1"] - df["p2"]
df["value_cond"] = pd.NA
df.loc[diff > 0, "value_cond"] = "left>right"
df.loc[diff < 0, "value_cond"] = "right>left"
df["value_cond3"] = "left=right"
df.loc[diff > 0, "value_cond3"] = "left>right"
df.loc[diff < 0, "value_cond3"] = "right>left"
df["value_cond"]  = pd.Categorical(df["value_cond"],  ["left>right","right>left"], ordered=True)
df["value_cond3"] = pd.Categorical(df["value_cond3"], ["left>right","left=right","right>left"], ordered=True)

# long table format ---------------------------------------------------------------------

dur_long = pd.concat([
    pd.DataFrame({
        "sub_id": df["sub_id"], "role": "First",
        "opt_type": df["FirstFixLoc"].map(loc_to_side),
        "dur_ms": df["FirstFixDur"],
        "value_cond": df["value_cond"],
        "OV_2": df["OV_2"],
    }),
    pd.DataFrame({
        "sub_id": df["sub_id"], "role": "Middle",
        "opt_type": df["MiddleDominantLoc"].map(loc_to_side),
        "dur_ms": df["eachMiddleFixDur"],
        "value_cond": df["value_cond"],
        "OV_2": df["OV_2"],
    }),
    pd.DataFrame({
        "sub_id": df["sub_id"], "role": "Final",
        "opt_type": df["FinalFixLoc"].map(loc_to_side),
        "dur_ms": df["FinalFixDur"],
        "value_cond": df["value_cond"],
        "OV_2": df["OV_2"],
    }),
], ignore_index=True)

dur_long = dur_long[
    dur_long["opt_type"].isin(["left","right"]) &
    (~dur_long["dur_ms"].isna()) &
    (dur_long["dur_ms"] > 0) &
    (~dur_long["value_cond"].isna()) &
    (~dur_long["OV_2"].isna())
].copy()

dur_long["role"]     = pd.Categorical(dur_long["role"], ["First","Middle","Final"], ordered=True)
dur_long["opt_type"] = pd.Categorical(dur_long["opt_type"], ["left","right"], ordered=True)
dur_long["log_dur"] = np.log(dur_long["dur_ms"])

# Locations (trial-level long)
loc_long = pd.concat([
    pd.DataFrame({
        "sub_id": df["sub_id"], "role": "First",
        "opt_type": df["FirstFixLoc"].map(loc_to_side),
        "value_cond3": df["value_cond3"],
        "OV_2": df["OV_2"],
    }),
    pd.DataFrame({
        "sub_id": df["sub_id"], "role": "Middle",
        "opt_type": df["MiddleDominantLoc"].map(loc_to_side),
        "value_cond3": df["value_cond3"],
        "OV_2": df["OV_2"],
    }),
    pd.DataFrame({
        "sub_id": df["sub_id"], "role": "Final",
        "opt_type": df["FinalFixLoc"].map(loc_to_side),
        "value_cond3": df["value_cond3"],
        "OV_2": df["OV_2"],
    }),
], ignore_index=True)

loc_long = loc_long[
    loc_long["opt_type"].isin(["left","right"]) &
    (~loc_long["value_cond3"].isna()) &
    (~loc_long["OV_2"].isna())
].copy()

loc_long["role"]     = pd.Categorical(loc_long["role"], ["First","Middle","Final"], ordered=True)
loc_long["opt_type"] = pd.Categorical(loc_long["opt_type"], ["left","right"], ordered=True)

# duration role*val cond* OV_2* opt_type ----------------
subj_ev_ov = (
    dur_long
    .groupby(["sub_id","role","value_cond","OV_2","opt_type"], observed=True, as_index=False)
    .agg(mean_ms=("dur_ms","mean"),
         median_ms=("dur_ms","median"),
         mean_log=("log_dur","mean"),
         sd_log=("log_dur","std"),
         n_trials=("dur_ms","size"))
)

dur_ev_ov_sum = (
    subj_ev_ov
    .groupby(["role","value_cond","OV_2","opt_type"], observed=True, as_index=False)
    .agg(n_subjects=("sub_id","nunique"),
         n_trials_total=("n_trials","sum"),
         mean_ms=("mean_ms","mean"),
         sd_ms=("mean_ms", sd1),
         median_ms=("median_ms","median"),

         # log summaries
         mu_log=("mean_log","mean"),
         sigma_log=("mean_log", sd1))
    .sort_values(["role","value_cond","OV_2","opt_type"])
)
dur_ev_ov_sum["se_ms"]  = dur_ev_ov_sum["sd_ms"]  / np.sqrt(dur_ev_ov_sum["n_subjects"].clip(lower=1))
dur_ev_ov_sum["se_log"] = dur_ev_ov_sum["sigma_log"] / np.sqrt(dur_ev_ov_sum["n_subjects"].clip(lower=1))
dur_ev_ov_sum["geom_mean_ms"] = np.exp(dur_ev_ov_sum["mu_log"])
dur_ev_ov_sum["geom_sd_multiplier"] = np.exp(dur_ev_ov_sum["sigma_log"])
dur_ev_ov_sum.to_csv(os.path.join(OUT_DIR, FN_DUR_EV_OV), index=False)
print("Saved:", os.path.abspath(os.path.join(OUT_DIR, FN_DUR_EV_OV)))

# probability role*val cond* OV_2* opt_type ----------------

counts = (
    loc_long
    .groupby(["sub_id","role","value_cond3","OV_2","opt_type"])
    .size()
    .unstack(fill_value=0)
    .reindex(columns=["left","right"], fill_value=0)
    .reset_index()
)

counts["total"] = counts["left"] + counts["right"]
counts = counts[counts["total"] > 0].copy()
counts["prob_left"]  = counts["left"]  / counts["total"]
counts["prob_right"] = counts["right"] / counts["total"]

subj_probs = counts.melt(
    id_vars=["sub_id","role","value_cond3","OV_2"],
    value_vars=["prob_left","prob_right"],
    var_name="opt_type",
    value_name="prob"
)
subj_probs["opt_type"] = subj_probs["opt_type"].map({"prob_left":"left","prob_right":"right"})
subj_probs["opt_type"] = pd.Categorical(subj_probs["opt_type"], ["left","right"], ordered=True)

prop_ev_ov = (
    subj_probs
    .groupby(["role","value_cond3","OV_2","opt_type"], observed=True, as_index=False)
    .agg(n_subjects=("sub_id","nunique"),
         mean_prob=("prob","mean"),
         sd_prob=("prob", sd1),
         median_prob=("prob","median"))
    .sort_values(["role","value_cond3","OV_2","opt_type"])
)

prop_ev_ov["se_prob"] = prop_ev_ov["sd_prob"] / np.sqrt(prop_ev_ov["n_subjects"].clip(lower=1))
prop_ev_ov["mean_percent"] = 100 * prop_ev_ov["mean_prob"]
prop_ev_ov["se_percent"]   = 100 * prop_ev_ov["se_prob"]
prop_ev_ov["median_percent"] = 100 * prop_ev_ov["median_prob"]
prop_ev_ov.to_csv(os.path.join(OUT_DIR, FN_PROP_EV_OV), index=False)
print("Saved:", os.path.abspath(os.path.join(OUT_DIR, FN_PROP_EV_OV)))
