# -*- coding: utf-8 -*-
# Probability that fixation is on E vs S by value dominance in one big plot.
# For each fixation type: 4 bars = (E>S,E), (E>S,S), (S>E,E), (S>E,S)
# Aggregation: mean over subjects (per-subject probabilities), SE across subjects
# Saves CSV and one PNG (percent scale).

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import ttest_rel, t

# --- Inputs ---
file = r"D:/Aberdeen_Uni_June24/cap/THESIS/Garcia_Analysis/data/data_sets/GarciaParticipants_Eye_Response_Feed_Allfix_addm_OV_Abs_CCT.csv"
excluded_subs = {1, 4, 5, 6, 14, 99}   # OV study excludes, Garcia: 1, 4, 5, 6, 14, 99

output_csv_summary = "FixLoc_Prob_ES_by_ExpectedValue_Garcia.csv"
output_csv_stats   = "FixLoc_Prob_ES_by_ExpectedValue_OV_stats_paired_ttests_Garcia.csv"
output_plot        = "fixation_location_probability_ES_by_value_Garcia.png"

# --- Helpers ---
def to_numeric(x):
    try:
        if pd.isna(x): return np.nan
        return float(str(x).strip())
    except Exception:
        return np.nan

def loc_to_opt(loc_val):
    if pd.isna(loc_val): return np.nan
    v = int(round(loc_val))
    return "E" if v == 1 else ("S" if v == 2 else np.nan)

def sd1(x): 
    x = np.asarray(x, dtype=float)
    return float(np.std(x, ddof=1)) if x.size > 1 else np.nan

# --- Load & checks ---
df = pd.read_csv(file)
required = ["sub_id","phase","p1","p2","FirstFixLoc","MiddleDominantLoc","FinalFixLoc"]
missing = [c for c in required if c not in df.columns]
if missing:
    raise ValueError(f"Missing columns: {missing}")

# --- Clean & filter ---
df["sub_id"] = pd.to_numeric(df["sub_id"], errors="coerce").astype("Int64")
df["phase"]  = df["phase"].astype(str).str.upper()
for col in ["p1","p2","FirstFixLoc","MiddleDominantLoc","FinalFixLoc"]:
    df[col] = df[col].map(to_numeric)

df = df[(df["phase"] == "ES") & (~df["sub_id"].isin(excluded_subs))].copy()

# Value dominance: E>S (p1>p2), S>E (p2>p1); ties dropped
cond = np.where(df["p1"] > df["p2"], "E>S",
        np.where(df["p2"] > df["p1"], "S>E", None)).astype(object)

df["value_cond"] = pd.Categorical(cond, ["E>S","S>E"], ordered=True)
df = df[~df["value_cond"].isna()].copy()

# --- Long table of fixation locations ---
first  = pd.DataFrame({"sub_id": df["sub_id"], "role":"First",  "value_cond": df["value_cond"],
                       "opt_type": df["FirstFixLoc"].map(loc_to_opt)})
middle = pd.DataFrame({"sub_id": df["sub_id"], "role":"Middle", "value_cond": df["value_cond"],
                       "opt_type": df["MiddleDominantLoc"].map(loc_to_opt)})
final  = pd.DataFrame({"sub_id": df["sub_id"], "role":"Final",  "value_cond": df["value_cond"],
                       "opt_type": df["FinalFixLoc"].map(loc_to_opt)})

loc_long = pd.concat([first, middle, final], ignore_index=True)
loc_long = loc_long[loc_long["opt_type"].isin(["E","S"])].copy()
loc_long["role"] = pd.Categorical(loc_long["role"], ["First","Middle","Final"], ordered=True)
loc_long["opt_type"] = pd.Categorical(loc_long["opt_type"], ["E","S"], ordered=True)

# --- Per-subject probabilities within (role × value_cond) ---
counts = (loc_long.groupby(["sub_id","role","value_cond","opt_type"])
          .size().unstack(fill_value=0).reindex(columns=["E","S"], fill_value=0)).reset_index()
counts["total"]  = counts["E"] + counts["S"]
counts = counts[counts["total"] > 0].copy()
counts["prob_E"] = counts["E"] / counts["total"]
counts["prob_S"] = counts["S"] / counts["total"]

# Reshape to long for descriptives across subjects
subj_probs = counts.melt(id_vars=["sub_id","role","value_cond"],
                         value_vars=["prob_E","prob_S"],
                         var_name="opt_type", value_name="prob")
subj_probs["opt_type"] = subj_probs["opt_type"].map({"prob_E":"E","prob_S":"S"})
subj_probs["opt_type"] = pd.Categorical(subj_probs["opt_type"], ["E","S"], ordered=True)

# --- Across-subject descriptive summary (mean of subject probs, SE across subjects) ---
summary = (subj_probs.groupby(["role","value_cond","opt_type"], observed=True, as_index=False)
           .agg(n_subjects=("sub_id","nunique"),
                mean_prob=("prob","mean"),
                sd_prob=("prob", sd1),
                median_prob=("prob","median"))
           .sort_values(["role","value_cond","opt_type"]))
summary["se_prob"]      = summary["sd_prob"] / np.sqrt(summary["n_subjects"].clip(lower=1))
summary["mean_percent"] = summary["mean_prob"] * 100.0
summary["se_percent"]   = summary["se_prob"]   * 100.0
summary["median_percent"] = summary["median_prob"] * 100.0

summary.to_csv(output_csv_summary, index=False)
print(f"Saved summary CSV → {os.path.abspath(output_csv_summary)}")

# --- Paired t-tests: E vs S probabilities within each role × value_cond (per-subject) ---
rows = []
for role in ["First","Middle","Final"]:
    for val_cond in ["E>S","S>E"]:
        subdat = counts[(counts["role"]==role) & (counts["value_cond"]==val_cond)].copy()
        if subdat.empty:
            continue

        # require paired E/S probabilities per subject
        # (counts already has E, S, total, prob_E, prob_S per subject in this bin)
        n = subdat.shape[0]
        if n < 2:
            continue

        mean_E = float(subdat["prob_E"].mean())
        mean_S = float(subdat["prob_S"].mean())

        diffs = subdat["prob_E"] - subdat["prob_S"]
        mean_diff = float(diffs.mean())
        sd_diff   = sd1(diffs)
        se_diff   = sd_diff / np.sqrt(n)

        # Paired t-test
        tval, pval = ttest_rel(subdat["prob_E"], subdat["prob_S"])
        dfree = n - 1

        # 95% CI for mean difference
        tcrit = t.ppf(0.975, dfree)
        ci_low  = mean_diff - tcrit * se_diff
        ci_high = mean_diff + tcrit * se_diff

        # Cohen's dz for paired data
        cohen_dz = mean_diff / sd_diff if sd_diff and not np.isnan(sd_diff) else np.nan

        # Optional context: total fix counts across subjects in this bin
        n_fix_E = int(subdat["E"].sum())
        n_fix_S = int(subdat["S"].sum())

        rows.append(dict(
            role=role, value_cond=val_cond,
            n_subjects=n,
            mean_prob_E=mean_E, mean_prob_S=mean_S,
            mean_diff_prob=mean_diff, sd_diff_prob=sd_diff, se_diff_prob=se_diff,
            t_value=float(tval), df=int(dfree), p_value=float(pval),
            ci_diff_low=ci_low, ci_diff_high=ci_high,
            cohen_dz=cohen_dz,
            total_fix_E=n_fix_E, total_fix_S=n_fix_S
        ))

stats_df = pd.DataFrame(rows).sort_values(["role","value_cond"])
stats_df.to_csv(output_csv_stats, index=False)
print(f"Saved stats CSV → {os.path.abspath(output_csv_stats)}")

# --- Plot: single figure (3 groups × 4 bars per group) ---
from matplotlib.patches import Patch

roles = ["First","Middle","Final"]
mapping = [("E>S","E"), ("E>S","S"), ("S>E","E"), ("S>E","S")]

mean_mat = np.full((len(roles), 4), np.nan)
se_mat   = np.full_like(mean_mat, np.nan, dtype=float)

for i, r in enumerate(roles):
    for k, (c,o) in enumerate(mapping):
        row = summary[(summary["role"]==r) & (summary["value_cond"]==c) & (summary["opt_type"]==o)]
        if len(row) == 1:
            mean_mat[i, k] = row["mean_percent"].values[0]
            se_mat[i, k]   = row["se_percent"].values[0]

x = np.arange(len(roles))
width = 0.18
pos = [x - 1.5*width, x - 0.5*width, x + 0.5*width, x + 1.5*width]

fig, ax = plt.subplots(figsize=(10, 5.8), dpi=150)

ax.bar(pos[0], mean_mat[:,0], width, yerr=np.nan_to_num(se_mat[:,0]), capsize=4, color="deepskyblue", edgecolor="black")
ax.bar(pos[1], mean_mat[:,1], width, yerr=np.nan_to_num(se_mat[:,1]), capsize=4, color="darkorchid", edgecolor="black")
ax.bar(pos[2], mean_mat[:,2], width, yerr=np.nan_to_num(se_mat[:,2]), capsize=4, color="deepskyblue", hatch="//", edgecolor="black")
ax.bar(pos[3], mean_mat[:,3], width, yerr=np.nan_to_num(se_mat[:,3]), capsize=4, color="darkorchid",  hatch="//", edgecolor="black")

ax.set_xticks(x, roles)
ax.set_ylabel("Probability (%)")
ax.set_ylim(0, 100)
ax.set_title("P(Fixate E vs S) by Expected Value and Option (mean ± SE)", pad=10)
ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)

handles = [
    Patch(facecolor="deepskyblue", edgecolor="black", label="E"),
    Patch(facecolor="darkorchid",  edgecolor="black", label="S"),
    Patch(facecolor="white", edgecolor="black", label="Solid = E>S"),
    Patch(facecolor="white", edgecolor="black", hatch="//", label="Hatched = S>E"),
]
ax.legend(handles=handles, loc="upper right", frameon=False,
          prop={"size": 9}, labelspacing=0.3, borderpad=0.3, handletextpad=0.6, borderaxespad=0.3)

fig.tight_layout()
fig.savefig(output_plot, dpi=300, bbox_inches="tight")
plt.show()

# --- Console summary for quick copy/paste ---
print("\nPaired t-tests (prob_E vs prob_S) within each role × value_cond:")
print(stats_df.to_string(index=False))