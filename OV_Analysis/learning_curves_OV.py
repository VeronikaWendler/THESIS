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

titles = {0: "Condition 90/70", 
          1: "Condition 80/40", 
          2: "Condition 60/20", 
          3: "Condition 30/10"}


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

g.figure.suptitle("Learning curves in the LE phase of EXP2 (mean ± 95 % CI)", y=1.02)
plt.tight_layout()
#plt.show()

plt.savefig("D:/Aberdeen_Uni_June24/cap/THESIS/OV_Analysis/Figures/learning_curves_LE_phase.png", dpi=600, bbox_inches='tight')

plt.show()

