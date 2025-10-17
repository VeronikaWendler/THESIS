# # Confidence ratings
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import ttest_rel
from statsmodels.stats.multitest import multipletests
import os


file_path = r'D:/Aberdeen_Uni_June24/cap/THESIS/Garcia_Analysis/data/data_sets/GarciaParticipants_Eye_Response_Feed_Allfix_addm_OV_Abs_CCT.csv'
data = pd.read_csv(file_path)

out_dir = r"D:/Aberdeen_Uni_June24/cap/THESIS/Garcia_Analysis/stats_TingGluth/confidence_ratings"


exclude_subjects = [6, 99]
data = data[(data['phase'] == 'SP') & (~data['sub_id'].isin(exclude_subjects))]

data['p1'] = pd.to_numeric(data['p1'], errors='coerce')
data['p1'] = data['p1'] / 100  # convert from 90.0 → 0.90

data['cho'] = pd.to_numeric(data['cho'], errors='coerce')

data = data.dropna(subset=['p1', 'cho', 'op1'])

probability_levels = np.arange(0.1, 1.0, 0.1)
tol = 0.01  


statistical_results = []

# paired t-tests for each probability level
for prob in probability_levels:
    participant_means_E, participant_means_S = [], []
    
    for participant in data['sub_id'].unique():
        subset_E = data[(np.abs(data['p1'] - prob) < tol) & (data['op1'] == 'E') & (data['sub_id'] == participant)]['cho']
        subset_S = data[(np.abs(data['p1'] - prob) < tol) & (data['op1'] == 'S') & (data['sub_id'] == participant)]['cho']
        
        if not subset_E.empty and not subset_S.empty:
            participant_means_E.append(subset_E.mean())
            participant_means_S.append(subset_S.mean())
    
    if len(participant_means_E) > 1 and len(participant_means_S) > 1:
        t_stat, p_val = ttest_rel(participant_means_E, participant_means_S)
        cohen_d = (np.mean(participant_means_E) - np.mean(participant_means_S)) / np.std(participant_means_E + participant_means_S, ddof=1)
        df = len(participant_means_E) - 1
        statistical_results.append((prob, df, t_stat, p_val, cohen_d))

stats_df = pd.DataFrame(statistical_results, columns=['Objective Probability', 'df', 'T-Statistic', 'P-Value', 'Cohen_d'])
stats_df['P-Value Corrected'] = multipletests(stats_df['P-Value'], method='bonferroni')[1]

custom_palette = {"E": "deepskyblue", "S": "darkorchid"}

plt.figure(figsize=(8, 6), facecolor='white')
sns.pointplot(data=data, x='p1', y='cho', hue='op1', dodge=True, capsize=0.1,
              markers=['o', 's'], linestyles=['-', '--'], palette=custom_palette)

plt.axhline(y=0.5, color='gray', linestyle='--', alpha=0.7)
plt.xlabel('Objective Probability', fontsize=18)
plt.ylabel('Rated Probability', fontsize=18)
plt.xticks(fontsize=14)
plt.yticks(fontsize=14)
plt.title('Stated Probabilities for E and S Options', fontsize=16)

legend = plt.legend(title="OP1", fontsize=16, title_fontsize=12)
legend.get_frame().set_facecolor('white')

plt.grid(False)
ax = plt.gca()
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
pdf_path = os.path.join(out_dir, 'stated_prob.png')
plt.savefig(pdf_path, format='png', dpi=300, bbox_inches='tight')
plt.show()

csv_path = os.path.join(out_dir, 'stats_results.csv')
stats_df.to_csv(csv_path, index=False)