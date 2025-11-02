import pandas as pd

# ---------------------- config ----------------------
file = r"D:/Aberdeen_Uni_June24/cap/THESIS/Garcia_Analysis/data/data_sets/GarciaParticipants_Eye_Response_Feed_Allfix_addm_OV_Abs_CCT.csv"
excluded_subs = {1, 4, 5, 6, 14, 99}

df = pd.read_csv(file)
df = df[~df['sub_id'].isin(excluded_subs)].copy()
df = df.dropna(subset=['OV_2', 'VD_2']).copy()

df['OV_2'] = df['OV_2'].astype(int)
df['VD_2'] = df['VD_2'].astype(int)

df['OV_VD_label'] = df['OV_2'].astype(str) + df['VD_2'].astype(str)
df['OV_VD_label'] = df['OV_VD_label'].astype(int)

combo_counts = df.groupby(['OV_2', 'VD_2']).size().reset_index(name='count')

print("\nOccurrences per OV_2, VD_2 combination ")
print(combo_counts)

print("\nValue counts of combined labels")
print(df['OV_VD_label'].value_counts().sort_index())

df.to_csv(file, index=False)
print(f"\n Saved updated file with OV_VD_label column to:\n{file}")
