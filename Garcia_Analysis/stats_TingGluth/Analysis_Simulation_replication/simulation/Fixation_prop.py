# Probability that fixation is on E vs S with SE bars (E blue, S red)

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

file = r"C:\Cluster_Github\HDDM_Vero\data_sets\data_sets_Garcia\GarciaParticipants_Eye_Response_Feed_Allfix_addm_OV_Abs_CCT.csv"
excluded_subs = {1, 4, 5, 6, 14, 99}
output_csv  = "Fix_Prop_ES.csv"   # same as before
output_plot = "fixation_location_probability_ES_bars_SE.png"       # NEW filename

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

df = pd.read_csv(file)
required = ["sub_id","phase","FirstFixLoc","FinalFixLoc","MiddleDominantLoc"]
missing = [c for c in required if c not in df.columns]
if missing:
    raise ValueError(f"Missing required columns: {missing}")

df["sub_id"] = pd.to_numeric(df["sub_id"], errors="coerce").astype("Int64")
df["phase"]  = df["phase"].astype(str).str.upper()
for col in ["FirstFixLoc","FinalFixLoc","MiddleDominantLoc"]:
    df[col] = df[col].map(to_numeric)

df = df[(df["phase"]=="ES") & (~df["sub_id"].isin(excluded_subs))].copy()

first  = pd.DataFrame({"sub_id": df["sub_id"], "role":"First",
                       "opt_type": df["FirstFixLoc"].map(loc_to_opt)})
middle = pd.DataFrame({"sub_id": df["sub_id"], "role":"Middle",
                       "opt_type": df["MiddleDominantLoc"].map(loc_to_opt)})
final  = pd.DataFrame({"sub_id": df["sub_id"], "role":"Final",
                       "opt_type": df["FinalFixLoc"].map(loc_to_opt)})

loc_long = pd.concat([first, middle, final], ignore_index=True)
loc_long = loc_long[loc_long["opt_type"].isin(["E","S"])].copy()
loc_long["role"] = pd.Categorical(loc_long["role"], ["First","Middle","Final"], ordered=True)
loc_long["opt_type"] = pd.Categorical(loc_long["opt_type"], ["E","S"], ordered=True)

# Per-subject probabilities by role
counts = (loc_long
          .groupby(["sub_id","role","opt_type"])
          .size()
          .unstack(fill_value=0)
          .reindex(columns=["E","S"], fill_value=0))

counts = counts.reset_index()
counts["total"] = counts["E"] + counts["S"]
counts = counts[counts["total"] > 0].copy()
counts["prob_E"] = counts["E"] / counts["total"]
counts["prob_S"] = counts["S"] / counts["total"]

subj_probs = counts.melt(
    id_vars=["sub_id","role"],
    value_vars=["prob_E","prob_S"],
    var_name="opt_type",
    value_name="prob"
)
subj_probs["opt_type"] = subj_probs["opt_type"].map({"prob_E":"E","prob_S":"S"})
subj_probs["opt_type"] = pd.Categorical(subj_probs["opt_type"], ["E","S"], ordered=True)

def sd1(x):
    return float(np.std(x, ddof=1)) if len(x)>1 else np.nan

summary = (subj_probs
    .groupby(["role","opt_type"], observed=True, as_index=False)
    .agg(n_subjects=("sub_id","nunique"),
         mean_prob=("prob","mean"),
         sd_prob=("prob", sd1),
         median_prob=("prob","median"))
    .sort_values(["role","opt_type"])
)

summary["se_prob"]      = summary["sd_prob"] / np.sqrt(summary["n_subjects"].clip(lower=1))
summary["mean_percent"] = summary["mean_prob"] * 100.0
summary["se_percent"]   = summary["se_prob"]   * 100.0
summary["median_percent"] = summary["median_prob"] * 100.0

summary.to_csv(output_csv, index=False)
print(f"Saved: {os.path.abspath(output_csv)}")

# Plot with SE bars 
roles     = ["First","Middle","Final"]
opt_types = ["E","S"]
role_idx  = {r:i for i,r in enumerate(roles)}

mean_grid = np.full((len(roles), len(opt_types)), np.nan)
se_grid   = np.full_like(mean_grid, np.nan, dtype=float)

for _, row in summary.iterrows():
    i = role_idx[row["role"]]
    j = opt_types.index(row["opt_type"])
    mean_grid[i, j] = row["mean_percent"]
    se_grid[i, j]   = row["se_percent"]

x = np.arange(len(roles))
width = 0.35

fig, ax = plt.subplots(figsize=(8,5), dpi=150)
ax.bar(x - width/2, mean_grid[:,0], width, label="E", color="deepskyblue",
       yerr=np.nan_to_num(se_grid[:,0]), capsize=4)
ax.bar(x + width/2, mean_grid[:,1], width, label="S", color="darkorchid",
       yerr=np.nan_to_num(se_grid[:,1]), capsize=4)

ax.set_ylabel("Probability (%)")
ax.set_xticks(x, roles)
ax.set_ylim(0, 100)
ax.set_title("Probability Fixation Falls on E vs S (mean ± SE)")
ax.legend(loc="upper center", ncols=2)
ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
fig.tight_layout()
fig.savefig(output_plot)
print(f"Saved: {os.path.abspath(output_plot)}")
plt.show()
