# Veronika Wendler
# computing learning curves


import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import sem                       

data = pd.read_csv("D:/Aberdeen_Uni_June24/cap/THESIS/OV_Analysis/data/data_sets/OVParticipants_Eye_Response_Feed_Allfix_addm_OV_Abs_CCT.csv")

data = data[data["phase"] == "LE"].copy()
data["trial_2"]  = data.groupby(["sub_id", "cond"]).cumcount() + 1
data["corr_pct"] = data["corr"] * 100          # seaborn wants raw units

titles = {0: "Condition 70/90", 
          1: "Condition 40/80", 
          2: "Condition 20/60", 
          3: "Condition 10/30"}

g = sns.relplot(
        data=data,
        x="trial_2",
        y="corr_pct",
        col="cond",
        col_wrap=2,
        kind="line",
        estimator="mean",
        ci=95,
        n_boot=1000,
        color="steelblue",           
        lw=2,
        height=3.5, 
        aspect=1.2,
        legend=False)               

for cond, ax in zip(sorted(data["cond"].unique()), g.axes.flat):
    ax.axhline(50, ls="--", lw=1, color="red")        
    ax.set_xticks([0, 10, 20, 30, 40])
    ax.set_xlim(0, 40)
    ax.set_ylim(0, 100)
    ax.set_xlabel("Trial")
    ax.set_ylabel("Accuracy (%)")
    ax.set_title(titles.get(cond, f"Cond {cond}"))    

g.figure.suptitle("Learning curves in the LE phase of EXP1 (mean ± 95 % CI)", y=1.02)
plt.tight_layout()
#plt.show()

plt.savefig("D:/Aberdeen_Uni_June24/cap/THESIS/Garcia_Analysis/Figures/learning_curves_LE_phase.png", dpi=600, bbox_inches='tight')

plt.show()







# # filter for LE phase
# data = data[data["phase"] == "LE"].copy()

# data["trial_2"] = data.groupby(["sub_id", "cond"]).cumcount() + 1

# assert data["trial_2"].max() == 40

# # Aggregate across participants
# curve = (data.groupby(["cond", "trial_2"])["corr"].agg(mean="mean", n="size").reset_index())

# # Standard error & 95 % CI for a proportion
# curve["sem"]  = np.sqrt(curve["mean"] * (1 - curve["mean"]) / curve["n"])
# curve["ci95"] = 1.96 * curve["sem"]

# # 5.  Plot – one panel per condition (0-3)
# fig, axes = plt.subplots(2, 2, figsize=(10, 8), sharex=True, sharey=True)
# axes = axes.ravel()

# for c, ax in enumerate(axes):
#     sub = curve[curve["cond"] == c]
#     ax.plot(sub["trial_2"], sub["mean"]*100, lw=2, label=f"Cond {c}")
#     ax.fill_between(sub["trial_2"],
#                     (sub["mean"]-sub["ci95"])*100,
#                     (sub["mean"]+sub["ci95"])*100,
#                     alpha=.25)
#     ax.axhline(50, ls="--", lw=1, color="red")        # chance level
#     ax.set_title(f"Condition {c}")
#     ax.set_xlim(1, 40)
#     ax.set_ylim(0, 100)
#     ax.set_xlabel("Trial")
#     ax.set_ylabel("Accuracy (%)")

# fig.suptitle("Learning curves for each condition (mean ± 95 % CI)")
# fig.tight_layout()
# plt.show()

# plt.savefig(
#     "D:/Aberdeen_Uni_June24/cap/THESIS/Garcia_Analysis/Figures/learning_curves_LE_phase.png",
#     dpi=600,               # high resolution
#     bbox_inches='tight',   # removes extra white space
#     transparent=False
# )