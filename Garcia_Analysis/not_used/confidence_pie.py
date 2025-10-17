import pandas as pd
import glob
import matplotlib.pyplot as plt
import seaborn as sns
from statsmodels.stats.anova import AnovaRM
import pingouin as pg

# load files
file_pattern = r'D:/Aberdeen_Uni_June24/cap/THESIS/Garcia_Analysis/data/sub-*/beh/EXP4_Garcia_participant_*.csv'
df = pd.concat((pd.read_csv(f, sep=',') for f in glob.glob(file_pattern)), ignore_index=True)
out_dir = r"D:/Aberdeen_Uni_June24/cap/THESIS/Garcia_Analysis/stats_TingGluth/confidence_ratings"

#filter out Image trials
df = df[df['selectedImageNamesArrayEXP'].str.contains('Pie', na=False)]

# get per-participant mean confidence
grp = (
    df
    .groupby(['SubID','p1'])['confidenceLevelsArrayEXP']
    .mean()
    .reset_index(name='mean_conf')
    .query("p1 != 0.5")
)

# Repeated-measures ANOVA
aov = AnovaRM(grp, depvar='mean_conf', subject='SubID', within=['p1']).fit()
print(aov)

# Pairwise within-subject t-tests (Bonferroni-corrected)
posthocs = pg.pairwise_ttests(
    data=grp,
    dv='mean_conf',
    within='p1',
    subject='SubID',
    padjust='bonf'
)
print(posthocs)

#plot
grp = (
    df
    .groupby(['SubID','p1'])['confidenceLevelsArrayEXP']
    .mean()
    .reset_index(name='mean_conf')
)

# drop p1 == 0.5 (pies only)
grp = grp[grp['p1'] != 0.5]

#plot violin of the distribution of participant means at each p1
plt.figure(figsize=(8,5))
sns.violinplot(
    x='p1',
    y='mean_conf',
    data=grp,
    color="#AA4E73",
    order=sorted(grp['p1'].unique())
)
plt.xlabel('Objective probability (p1)')
plt.ylabel('Mean confidence level')
plt.title('Mean confidence for S-options (Pie Charts) by objective probability - EXP1')
plt.tight_layout()
plot_path = os.path.join(out_dir, 'mean_conf_violin_pie.png')
plt.savefig(plot_path, dpi=300)
plt.close()

# ANOVA text summary
anova_txt = os.path.join(out_dir, 'anova_results_pie.txt')
with open(anova_txt, 'w') as f:
    f.write(aov.summary().as_text())
print(f"Saved ANOVA summary to {anova_txt}")

# post-hoc to CSV
posthoc_csv = os.path.join(out_dir, 'posthoc_results_pie.csv')
posthocs.to_csv(posthoc_csv, index=False)
print(f"Saved post-hoc results to {posthoc_csv}")

# OPTIONAL: combine both into one Excel workbook
excel_path = os.path.join(out_dir, 'stats_results_pie.xlsx')
with pd.ExcelWriter(excel_path) as writer:
    # if your aov has anova_table attribute (Statsmodels >=0.14)
    try:
        aov.anova_table.to_excel(writer, sheet_name='ANOVA')
    except AttributeError:
        # fallback: write summary text into its own sheet
        df = pd.DataFrame([aov.summary().tables[0].data])
        df.to_excel(writer, sheet_name='ANOVA', index=False, header=False)
    posthocs.to_excel(writer, sheet_name='PostHoc', index=False)
print(f"Saved combined stats to {excel_path}")