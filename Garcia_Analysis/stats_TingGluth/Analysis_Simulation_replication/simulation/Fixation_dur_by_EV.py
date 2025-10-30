# Fixation durations (ms) by value
# For each fixation type (First/Middle/Final): 4 bars
# Aggregation: mean over subjects (mean of per-subject means), SE across subjects
# Saves CSV and one PNG

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import ttest_rel, t
from scipy.stats import t


file = r"D:/Aberdeen_Uni_June24/cap/THESIS/Garcia_Analysis/data/data_sets/GarciaParticipants_Eye_Response_Feed_Allfix_addm_OV_Abs_CCT.csv"
excluded_subs = {1,4,5,6,14,99}   # OV: 6, 14, 20, 26, 2, 9, 18   # Garcia: 1,4,5,6,14,99

output_csv_summary = "fixation_duration_ES_by_value_condition_subjectmean_Garcia.csv"
output_csv_stats   = "fixation_duration_ES_by_value_stats_paired_ttests_Garcia.csv"
output_plot        = "fixation_duration_ES_by_value_Garcia.png"

# Helpers 
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

# Load & checks 
df = pd.read_csv(file)
required = ["sub_id","phase","p1","p2",
            "FirstFixLoc","FirstFixDur",
            "MiddleDominantLoc","MiddleDominantDur",
            "FinalFixLoc","FinalFixDur"]
missing = [c for c in required if c not in df.columns]
if missing:
    raise ValueError(f"Missing columns: {missing}")

# Clean & filter 
df["sub_id"] = pd.to_numeric(df["sub_id"], errors="coerce").astype("Int64")
df["phase"]  = df["phase"].astype(str).str.upper()
for col in ["p1","p2","FirstFixLoc","MiddleDominantLoc","FinalFixLoc",
            "FirstFixDur","MiddleDominantDur","FinalFixDur"]:
    df[col] = df[col].map(to_numeric)

df = df[(df["phase"] == "ES") & (~df["sub_id"].isin(excluded_subs))].copy()

# Value dominance: E>S (p1>p2), S>E (p2>p1)
cond = np.where(df["p1"] > df["p2"], "E>S",
        np.where(df["p2"] > df["p1"], "S>E", None)).astype(object)

df["value_cond"] = pd.Categorical(cond, ["E>S","S>E"], ordered=True)
df = df[~df["value_cond"].isna()].copy()

# Long table: role × option × value_cond with durations 
first  = pd.DataFrame({"sub_id": df["sub_id"], "role":"First",  "value_cond": df["value_cond"],
                       "opt_type": df["FirstFixLoc"].map(loc_to_opt),        "dur_ms": df["FirstFixDur"]})
middle = pd.DataFrame({"sub_id": df["sub_id"], "role":"Middle", "value_cond": df["value_cond"],
                       "opt_type": df["MiddleDominantLoc"].map(loc_to_opt), "dur_ms": df["MiddleDominantDur"]})
final  = pd.DataFrame({"sub_id": df["sub_id"], "role":"Final",  "value_cond": df["value_cond"],
                       "opt_type": df["FinalFixLoc"].map(loc_to_opt),        "dur_ms": df["FinalFixDur"]})

dur_long = pd.concat([first, middle, final], ignore_index=True)
dur_long = dur_long[dur_long["opt_type"].isin(["E","S"]) & (~dur_long["dur_ms"].isna())].copy()
dur_long["role"] = pd.Categorical(dur_long["role"], ["First","Middle","Final"], ordered=True)
dur_long["opt_type"] = pd.Categorical(dur_long["opt_type"], ["E","S"], ordered=True)

# subject means
subj = (dur_long.groupby(["sub_id","role","value_cond","opt_type"], observed=True, as_index=False)
        .agg(mean_ms=("dur_ms","mean"),
             median_ms=("dur_ms","median"),
             n_trials=("dur_ms","size")))

# Across-subject descriptive summary
summary = (subj.groupby(["role","value_cond","opt_type"], observed=True, as_index=False)
           .agg(n_subjects=("sub_id","nunique"),
                n_trials_total=("n_trials","sum"),
                mean_ms=("mean_ms","mean"),
                sd_ms=("mean_ms", sd1),
                median_ms=("median_ms","median"))
           .sort_values(["role","value_cond","opt_type"]))
summary["se_ms"] = summary["sd_ms"] / np.sqrt(summary["n_subjects"].clip(lower=1))

summary.to_csv(output_csv_summary, index=False)
print(f"Saved summary CSV → {os.path.abspath(output_csv_summary)}")

# Paired t-tests: E vs S within each role × value_cond (per-subject means) 
rows = []
for role in ["First","Middle","Final"]:
    for val_cond in ["E>S","S>E"]:
        subdat = subj[(subj["role"]==role) & (subj["value_cond"]==val_cond)].copy()
        if subdat.empty:
            continue
        wide = subdat.pivot(index="sub_id", columns="opt_type", values="mean_ms")

        if not set(["E","S"]).issubset(wide.columns):
            continue
        wide = wide.dropna(subset=["E","S"])
        n = wide.shape[0]
        if n < 2:
            continue

        mean_E = float(wide["E"].mean())
        mean_S = float(wide["S"].mean())

        diffs = wide["E"] - wide["S"]
        mean_diff = float(diffs.mean())
        sd_diff   = sd1(diffs)
        se_diff   = sd_diff / np.sqrt(n)

        # Paired t-test
        tval, pval = ttest_rel(wide["E"], wide["S"])
        dfree = n - 1

        tcrit = t.ppf(0.975, dfree)
        ci_low  = mean_diff - tcrit * se_diff
        ci_high = mean_diff + tcrit * se_diff
        cohen_dz = mean_diff / sd_diff if sd_diff and not np.isnan(sd_diff) else np.nan

        n_trials_E = int(subdat.loc[subdat["opt_type"]=="E","n_trials"].sum())
        n_trials_S = int(subdat.loc[subdat["opt_type"]=="S","n_trials"].sum())

        rows.append(dict(
            role=role, value_cond=val_cond,
            n_subjects=n,
            mean_E_ms=mean_E, mean_S_ms=mean_S,
            mean_diff_ms=mean_diff, sd_diff_ms=sd_diff, se_diff_ms=se_diff,
            t_value=float(tval), df=int(dfree), p_value=float(pval),
            ci_diff_low_ms=ci_low, ci_diff_high_ms=ci_high,
            cohen_dz=cohen_dz,
            total_trials_E=n_trials_E, total_trials_S=n_trials_S
        ))

stats_df = pd.DataFrame(rows).sort_values(["role","value_cond"])
stats_df.to_csv(output_csv_stats, index=False)
print(f"Saved stats CSV → {os.path.abspath(output_csv_stats)}")

from matplotlib.patches import Patch

roles = ["First","Middle","Final"]
mapping = [("E>S","E"), ("E>S","S"), ("S>E","E"), ("S>E","S")]

mean_mat = np.full((len(roles), 4), np.nan)
se_mat   = np.full_like(mean_mat, np.nan, dtype=float)

for i, r in enumerate(roles):
    for k, (c,o) in enumerate(mapping):
        row = summary[(summary["role"]==r) & (summary["value_cond"]==c) & (summary["opt_type"]==o)]
        if len(row) == 1:
            mean_mat[i, k] = row["mean_ms"].values[0]
            se_mat[i, k]   = row["se_ms"].values[0]

x = np.arange(len(roles))
width = 0.18
pos = [x - 1.5*width, x - 0.5*width, x + 0.5*width, x + 1.5*width]

fig, ax = plt.subplots(figsize=(10, 5.8), dpi=150)

ax.bar(pos[0], mean_mat[:,0], width, yerr=np.nan_to_num(se_mat[:,0]), capsize=4, color="deepskyblue", edgecolor="black")
ax.bar(pos[1], mean_mat[:,1], width, yerr=np.nan_to_num(se_mat[:,1]), capsize=4, color="darkorchid", edgecolor="black")
ax.bar(pos[2], mean_mat[:,2], width, yerr=np.nan_to_num(se_mat[:,2]), capsize=4, color="deepskyblue", hatch="//", edgecolor="black")
ax.bar(pos[3], mean_mat[:,3], width, yerr=np.nan_to_num(se_mat[:,3]), capsize=4, color="darkorchid",  hatch="//", edgecolor="black")

ax.set_xticks(x, roles)
ax.set_ylabel("Fixation Duration (ms)")
ax.set_title("Fixation Duration by Value Dominance and Option (mean ± SE)", pad=10)
ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)

handles = [
    Patch(facecolor="deepskyblue", edgecolor="black", label="E (option)"),
    Patch(facecolor="darkorchid",  edgecolor="black", label="S (option)"),
    Patch(facecolor="white", edgecolor="black", label="Solid = E>S"),
    Patch(facecolor="white", edgecolor="black", hatch="//", label="Hatched = S>E"),
]
ax.legend(handles=handles, loc="upper left", frameon=False,
          prop={"size": 9}, labelspacing=0.3, borderpad=0.3, handletextpad=0.6, borderaxespad=0.3)

fig.tight_layout()
fig.savefig(output_plot, dpi=300, bbox_inches="tight")
plt.show()

print("\nPaired t-tests (E vs S) within each role × value_cond:")
print(stats_df.to_string(index=False))


##########


subj_role = (dur_long.groupby(["sub_id", "role"], observed=True, as_index=False)
             .agg(mean_ms=("dur_ms", "mean"),
                  median_ms=("dur_ms", "median"),
                  n_trials=("dur_ms", "size")))

overall_role = (subj_role.groupby("role", observed=True, as_index=False)
                .agg(n_subjects=("sub_id", "nunique"),
                     n_trials_total=("n_trials", "sum"),
                     mean_ms=("mean_ms", "mean"),
                     sd_ms=("mean_ms", sd1),
                     median_ms=("median_ms", "median"))
                .sort_values("role"))

overall_role["se_ms"] = overall_role["sd_ms"] / np.sqrt(overall_role["n_subjects"].clip(lower=1))

def ci95(row):
    n = int(row["n_subjects"])
    if n <= 1 or np.isnan(row["sd_ms"]):
        return (np.nan, np.nan)
    tcrit = t.ppf(0.975, n-1)
    se = row["sd_ms"] / np.sqrt(n)
    return (row["mean_ms"] - tcrit * se, row["mean_ms"] + tcrit * se)

cis = overall_role.apply(ci95, axis=1, result_type="expand")
overall_role["ci_low_ms"]  = cis[0]
overall_role["ci_high_ms"] = cis[1]

# Save CSV
output_csv_overall = "fixation_duration_ES_overall_by_role_Garcia.csv"
overall_role.to_csv(output_csv_overall, index=False)
print(f"Saved {os.path.abspath(output_csv_overall)}")

# Console print for quick reference
print("\nOverall fixation durations by role:")
print(overall_role[["role","n_subjects","mean_ms","sd_ms","se_ms","ci_low_ms","ci_high_ms","n_trials_total"]]
      .to_string(index=False))