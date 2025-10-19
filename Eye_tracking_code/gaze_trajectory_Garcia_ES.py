# Veronika Wendler
# 19.06.2025
# program that creates gaze trajectories during the learning phase (similar to plots in Checchi et al.,2025) to understand dynamic gaze allocation from stimuöi preseontation to outcome fixation

#libraries
import pandas as pd
import numpy as np 
import ast
import matplotlib.pyplot as plt
import itertools
from matplotlib.lines import Line2D
import os
from collections import defaultdict
from scipy.stats import ttest_rel
from scipy.stats import ttest_rel
from statsmodels.stats.multitest import multipletests
import ast 



# big file that has everything
df_raw = pd.read_csv("D:/Aberdeen_Uni_June24/cap/THESIS/Garcia_Analysis/data/data_sets/GarciaParticipants_Eye_Response_Feed_Allfix_addm_OV_Abs_CCT.csv")

# figures dir
fig_dir = "D:/Aberdeen_Uni_June24/cap/THESIS/Garcia_Analysis/Figures/gaze_trajectory"
if not fig_dir:
    os.mkdir(fig_dir)

BIN_MS    = 25                       # bin width in ms
ROI   = {1: 'left', 2: 'right'}      # ROI

raw = df_raw.query("phase == 'ES'").copy()   # keep only Learning‐phase trials
# Exclude bad subjects globally
EXCLUDE_SUBS = {1, 4, 5, 6, 14, 99}
raw = raw[~raw["sub_id"].isin(EXCLUDE_SUBS)].copy()



rt_glob_mean = raw.rtime.mean() * 1000          
rt_glob_sem  = raw.rtime.std(ddof=1) / np.sqrt(len(raw)) * 1000
raw['corr'] = raw['corr'].astype(int)


# helper for response phase -----------------------------------------------------------------------

def trial_bins(row):
    try:
        fixes = ast.literal_eval(row.Fixations) or []
    except Exception:
        return None
    if not fixes:
        return None

    hi_side = 'left' if row.p1 > row.p2 else 'right'
    t0      = fixes[0]['start_time']
    bins    = {}

    for f in fixes:
        side   = ROI.get(f['roi'],'none')
        is_hi  = side == hi_side
        for b in range(f['start_time']//BIN_MS - t0//BIN_MS,
                       (f['end_time']-1)//BIN_MS - t0//BIN_MS + 1):
            if b not in bins:
                bins[b] = [0,0]          # [high, low]
            bins[b][0 if is_hi else 1] += 1

    out = pd.DataFrame.from_dict(bins, orient='index', columns=['hi','lo'])
    out['sub'] = row.sub_id
    out['trial'] = row.TrialID
    out['bin'] = out.index
    out['hi_ratio'] = out.hi / (out.hi+out.lo)
    return out[['sub','trial','bin','hi_ratio']]


# helper for ES phase -----------------------------------------------------------------------
def trial_bins_ES(row):
    try:
        fixes = ast.literal_eval(row.Fixations) or []
    except Exception:
        return None
    if not fixes:
        return None

    E_side = 'left'
    t0 = fixes[0]['start_time']
    bins = {}

    for f in fixes:
        side  = ROI.get(f['roi'], 'none')
        is_E = side == E_side
        for b in range(f['start_time']//BIN_MS - t0//BIN_MS,
                       (f['end_time']-1)//BIN_MS - t0//BIN_MS + 1):
            if b not in bins:
                bins[b] = [0, 0]
            bins[b][0 if is_E else 1] += 1

    out = pd.DataFrame.from_dict(bins, orient='index', columns=['E','S'])
    out['sub'] = row.sub_id
    out['trial'] = row.TrialID
    out['bin'] = out.index
    out['ES_ratio'] = out.E / (out.E+out.S)
    return out[['sub','trial','bin','ES_ratio']]


def trial_bins_EvS_byChoice(row):
    #Return per-25 ms bins:
    #E_cnt, S_cnt
    # identify trial group where E-chosen (cho==1) or S-chosen (cho==2)
    
    try:
        fixes = ast.literal_eval(row.Fixations) or []
    except Exception:
        return None
    if not fixes:
        return None

    E_side, S_side = 'left', 'right'
    t0, bins = fixes[0]['start_time'], {}

    for f in fixes:
        side = ROI.get(f['roi'], 'none')            
        idx  = 0 if side == E_side else 1           
        for b in range(f['start_time']//BIN_MS - t0//BIN_MS,
                       (f['end_time']-1)//BIN_MS - t0//BIN_MS + 1):
            if b not in bins:
                bins[b] = [0, 0]                   
            bins[b][idx] += 1

    out = pd.DataFrame.from_dict(bins, orient='index',
                                 columns=['E_cnt', 'S_cnt'])
    out['sub']      = row.sub_id
    out['trial']    = row.TrialID
    out['bin']      = out.index
    out['cho_grp']  = 'E-chosen' if row.cho == 1 else 'S-chosen'
    out['ratio_E']  = out.E_cnt / (out.E_cnt + out.S_cnt)
    out['ratio_S']  = 1 - out['ratio_E']
    return out[['sub','trial','bin','cho_grp','ratio_E','ratio_S']]




def trial_bins_EvS_byAcc(row):
    #For every 25-ms bin during ES trials count E- and S-fixations.
    #Attach a trial label:
    # Correct    (row.corr == 1)
    # Incorrect  (row.corr == 0)
    
    try:
        fixes = ast.literal_eval(row.Fixations) or []
    except Exception:
        return None
    if not fixes:
        return None

    E_side, S_side = 'left', 'right'
    t0, bins = fixes[0]['start_time'], {}

    for f in fixes:
        side = ROI.get(f['roi'], 'none')          # 'left' / 'right' / 'none'
        idx  = 0 if side == E_side else 1         # 0 = E, 1 = S
        for b in range(f['start_time']//BIN_MS - t0//BIN_MS,
                       (f['end_time']-1)//BIN_MS - t0//BIN_MS + 1):
            if b not in bins:
                bins[b] = [0, 0]                  # [E_cnt, S_cnt]
            bins[b][idx] += 1

    out = pd.DataFrame.from_dict(bins, orient='index',
                                 columns=['E_cnt', 'S_cnt'])
    out['sub']   = row.sub_id
    out['trial'] = row.TrialID
    out['bin']   = out.index
    out['acc_grp'] = 'Correct' if row['corr'] == 1 else 'Incorrect'
    out['ratio_E'] = out.E_cnt / (out.E_cnt + out.S_cnt)
    out['ratio_S'] = 1 - out['ratio_E']
    return out[['sub','trial','bin','acc_grp','ratio_E','ratio_S']]



# bin table per trial 
tbl = raw.apply(trial_bins, axis=1).dropna()
tbl = pd.concat(tbl.tolist(), ignore_index=True)
rts = (raw[['sub_id','TrialID','rtime']]
       .rename(columns={'sub_id':'sub','TrialID':'trial','rtime':'rt'}))
rts['rt_bin'] = (rts.rt*1000 // BIN_MS).astype(int)
tbl = tbl.merge(rts, on=['sub','trial'])
tbl = tbl[tbl.bin <= tbl.rt_bin]         

# time in ms relative to onset
tbl['time_ms'] = tbl['bin'] * BIN_MS
# subject-average trajectory 
sub_traj = (tbl.groupby(['sub','time_ms']).hi_ratio.mean().reset_index(name='hi_ratio'))

# grand mean & SEM across participants
grand = (sub_traj.groupby('time_ms')
                   .agg(mean=('hi_ratio','mean'),
                        sem =('hi_ratio', lambda x: x.std()/np.sqrt(len(x))))
                   .reset_index())
grand['lo_mean'] = 1 - grand['mean']   
grand['lo_sem' ] = grand['sem']        

rt_mean = raw.rtime.mean()*1000            
rt_sem  = raw.rtime.std(ddof=1)/np.sqrt(len(raw))*1000

# table for ES phase
ES_tbl = raw.apply(trial_bins_ES, axis=1).dropna()
ES_tbl = pd.concat(ES_tbl.tolist(), ignore_index=True)
rts_ES = (raw[['sub_id','TrialID','rtime']]
       .rename(columns={'sub_id':'sub','TrialID':'trial','rtime':'rt'}))
rts_ES['rtES_bin'] = (rts_ES.rt*1000 // BIN_MS).astype(int)
ES_tbl = ES_tbl.merge(rts_ES, on=['sub','trial'])
ES_tbl = ES_tbl[ES_tbl.bin <= ES_tbl.rtES_bin]           
ES_tbl['timeES_ms'] = ES_tbl['bin'] * BIN_MS


# time windows analysis
windows = [(0, 300), (300, 600), (600, 900), (900, 1200)]
results = []

for (w_start, w_end) in windows:
    # average within window per subject
    df_win = (ES_tbl.query(f"{w_start} <= timeES_ms < {w_end}")
                      .groupby("sub")
                      .ES_ratio.mean()
                      .reset_index())
    df_win["S_ratio"] = 1 - df_win["ES_ratio"]

    # paired t-test across subjects
    tval, pval = ttest_rel(df_win["ES_ratio"], df_win["S_ratio"])

    # SEM across subjects
    E_sem = df_win["ES_ratio"].std(ddof=1) / np.sqrt(len(df_win))
    S_sem = df_win["S_ratio"].std(ddof=1) / np.sqrt(len(df_win))

    results.append({
        "window": f"{w_start}-{w_end} ms",
        "E_mean": df_win["ES_ratio"].mean(),
        "E_sem":  E_sem,
        "S_mean": df_win["S_ratio"].mean(),
        "S_sem":  S_sem,
        "diff": df_win["ES_ratio"].mean() - df_win["S_ratio"].mean(),
        "t": tval,
        "p": pval,
        "n_subs": df_win["sub"].nunique()
    })

results_df = pd.DataFrame(results)
print(results_df)
stats_path = os.path.join(fig_dir, "ES_EvS_window_stats.csv") 
results_df.to_csv(stats_path, index=False)



# average trajectory per participant
sub_trajES = (ES_tbl.groupby(['sub','timeES_ms']).ES_ratio.mean().reset_index(name='ES_ratio'))

# grand mean & SEM across participants
grandES = (sub_trajES.groupby('timeES_ms')
                   .agg(mean=('ES_ratio','mean'),
                        sem =('ES_ratio', lambda x: x.std()/np.sqrt(len(x))))
                   .reset_index())
grandES['S_mean'] = 1 - grandES['mean']  
grandES['S_sem' ] = grandES['sem']     
rtES_mean = raw.rtime.mean()*1000           
rtES_sem  = raw.rtime.std(ddof=1)/np.sqrt(len(raw))*1000

#----------------------------------------------------------------------------
# chocie E vs S

t_cap = int(rt_glob_mean + 100)    
choice_tbl = raw.apply(trial_bins_EvS_byChoice, axis=1).dropna()
choice_tbl = pd.concat(choice_tbl.tolist(), ignore_index=True)
choice_tbl['time_ms'] = choice_tbl['bin'] * BIN_MS
choice_tbl = choice_tbl[choice_tbl.time_ms <= rt_glob_mean + 100] 

sub_choice = (choice_tbl
              .groupby(['sub','cho_grp','time_ms'])
              .agg(E_mean=('ratio_E','mean'),
                   S_mean=('ratio_S','mean'))
              .reset_index())

grp_choice = (sub_choice
              .groupby(['cho_grp','time_ms'])
              .agg(E_m   = ('E_mean','mean'),
                   E_sem = ('E_mean', lambda x: x.std(ddof=1)/np.sqrt(len(x))),
                   S_m   = ('S_mean','mean'),
                   S_sem = ('S_mean', lambda x: x.std(ddof=1)/np.sqrt(len(x))))
              .reset_index())

# accuracy-based E vs S plot 
acc_tbl = raw.apply(trial_bins_EvS_byAcc, axis=1).dropna()
acc_tbl = pd.concat(acc_tbl.tolist(), ignore_index=True)
acc_tbl['time_ms'] = acc_tbl['bin'] * BIN_MS
acc_tbl = acc_tbl[acc_tbl.time_ms <= rt_glob_mean + 100] 

sub_acc = (acc_tbl
           .groupby(['sub','acc_grp','time_ms'])
           .agg(E_mean=('ratio_E','mean'),
                S_mean=('ratio_S','mean'))
           .reset_index())

grp_acc = (sub_acc
           .groupby(['acc_grp','time_ms'])
           .agg(E_m   = ('E_mean','mean'),
                E_sem = ('E_mean', lambda x: x.std(ddof=1)/np.sqrt(len(x))),
                S_m   = ('S_mean','mean'),
                S_sem = ('S_mean', lambda x: x.std(ddof=1)/np.sqrt(len(x))))
           .reset_index())

# helper for smaller legend
legend_kw = dict(loc='upper right', frameon=False,
                 fontsize=8,  title_fontsize=9,
                 borderpad=0.2, labelspacing=0.25, handlelength=2)
line_handles = [
    Line2D([], [], color='darkorchid',  lw=2, label='higher'),
    Line2D([], [], color='darkorchid', ls='--', lw=2, label='lower')
]
line_handles2 = [
    Line2D([], [], color='deepskyblue',  lw=2, label='E'),
    Line2D([], [], color='darkorchid', ls='--', lw=2, label='S')
]

# Choice figure - stimulus-locked ------------------------------------------

fig_choice, ax_c = plt.subplots(figsize=(4.5, 3))

# RT + mean SD
ax_c.axvspan(rt_mean-rt_sem, rt_mean+rt_sem,
             color='lightgrey', alpha=.5, zorder=0)
# RT mean line on top of the band
ax_c.axvline(rt_mean, color='grey', lw=1.2, zorder=1)

# higher-EV
ax_c.plot(grand.time_ms, grand['mean'],
          color='darkorchid', lw=2)
ax_c.fill_between(grand.time_ms,
                  grand['mean']-grand['sem'],
                  grand['mean']+grand['sem'],
                  color='darkorchid', alpha=.25, label='_nolegend_')

# lower-EV
ax_c.plot(grand.time_ms, grand['lo_mean'],
          color='darkorchid', ls='--', lw=2)
ax_c.fill_between(grand.time_ms,
                  grand['lo_mean']-grand['lo_sem'],
                  grand['lo_mean']+grand['lo_sem'],
                  color='darkorchid', alpha=.20, label='_nolegend_')

ax_c.axhline(.5, ls='--', color='k')
ax_c.axvline(0,  ls='--', color='k')
ax_c.set(xlim=(0, rt_mean+100), ylim=(0, 1),
         xlabel='Time (ms)', ylabel='Fixation ratio',
         title='Choice (stimulus-locked)')

ax_c.legend(handles=line_handles, title='Gaze on', **legend_kw)
fig_choice.tight_layout()
plot_path= os.path.join(fig_dir, "gaze_choice_ES_EV.png")
fig_choice.savefig(plot_path, dpi=300, bbox_inches='tight')
fig_choice.show()

# Feedback figure - stimulus-locked ------------------------------------------

fig_choiceES, ax_ES = plt.subplots(figsize=(4.5, 3))

# RT band - mean SD
ax_ES.axvspan(rtES_mean-rtES_sem, rtES_mean+rtES_sem,
             color='lightgrey', alpha=.5, zorder=0)
ax_ES.axvline(rtES_mean, color='grey', lw=1.2, zorder=1)

# higher-EV
ax_ES.plot(grandES.timeES_ms, grandES['mean'],
          color='deepskyblue', lw=2)
ax_ES.fill_between(grandES.timeES_ms,
                  grandES['mean']-grandES['sem'],
                  grandES['mean']+grandES['sem'],
                  color='deepskyblue', alpha=.25, label='_nolegend_')

# lower-EV
ax_ES.plot(grandES.timeES_ms, grandES['S_mean'],
          color='darkorchid', ls='--', lw=2)
ax_ES.fill_between(grandES.timeES_ms,
                  grandES['S_mean']-grandES['S_sem'],
                  grandES['S_mean']+grandES['S_sem'],
                  color='darkorchid', alpha=.20, label='_nolegend_')

ax_ES.axhline(.5, ls='--', color='k')
ax_ES.axvline(0,  ls='--', color='k')
ax_ES.set(xlim=(0, rtES_mean+100), ylim=(0, 1),
         xlabel='Time (ms)', ylabel='Fixation ratio',
         title='Choice E vs S (stimulus-locked)')
ax_ES.set_xlabel("Time (ms)", fontsize=13)
ax_ES.set_ylabel("Fixation ratio", fontsize=13)

ax_ES.legend(handles=line_handles2, title='Gaze on', **legend_kw)
fig_choiceES.tight_layout()
plot_path2= os.path.join(fig_dir, "gaze_choice_ES_ES.png")
fig_choiceES.savefig(plot_path2, dpi=300, bbox_inches='tight')
fig_choiceES.show()

# colours 
fig, ax = plt.subplots(figsize=(4.8,3.3))

# RT cue
ax.axvspan(rt_glob_mean-rt_glob_sem, rt_glob_mean+rt_glob_sem,
           color='lightgrey', alpha=.35, zorder=0)
ax.axvline(rt_glob_mean, color='grey', lw=1.2, zorder=1)

colour_grp = {'E-chosen':'deepskyblue', 'S-chosen':'darkorchid'}

for grp, df in grp_choice.groupby('cho_grp'):
    col = colour_grp[grp]

    # solid  = gaze to E
    ax.plot(df.time_ms, df.E_m,  color=col, lw=2,
            label=f'{grp} – gaze E')
    ax.fill_between(df.time_ms,
                    df.E_m-df.E_sem, df.E_m+df.E_sem,
                    color=col, alpha=.20)

    # dashed = gaze to S
    ax.plot(df.time_ms, df.S_m,  color=col, ls='--', lw=2,
            label=f'{grp} – gaze S')
    ax.fill_between(df.time_ms,
                    df.S_m-df.S_sem, df.S_m+df.S_sem,
                    color=col, alpha=.20)

ax.axhline(.5, ls='--', color='k')
ax.axvline(0,  ls='--', color='k')
ax.set(xlim=(0, rt_glob_mean+100), ylim=(0,1),
       xlabel='Time (ms)', ylabel='Fixation ratio',
       title='ES phase – gaze to E vs S, split by chosen item')
ax.legend(frameon=False, fontsize=7, ncol=1,
          bbox_to_anchor=(1.02,1), loc='upper left')
fig.tight_layout()
fig.savefig(os.path.join(fig_dir, 'ES_EvS_byChoice.png'),
            dpi=300, bbox_inches='tight')
plt.show()

#--------------------------------------

fig, ax = plt.subplots(figsize=(4.8,3.3))

# grey RT window + mean
ax.axvspan(rt_glob_mean-rt_glob_sem, rt_glob_mean+rt_glob_sem,
           color='lightgrey', alpha=.35, zorder=0)
ax.axvline(rt_glob_mean, color='grey', lw=1.2, zorder=1)

colour_grp = {'Correct':'seagreen', 'Incorrect':'crimson'}

for grp, df in grp_acc.groupby('acc_grp'):
    col = colour_grp[grp]

    # solid  = gaze to E
    ax.plot(df.time_ms, df.E_m,  color=col, lw=2,
            label=f'{grp} – gaze E')
    ax.fill_between(df.time_ms,
                    df.E_m-df.E_sem, df.E_m+df.E_sem,
                    color=col, alpha=.20)

    # dashed = gaze to S
    ax.plot(df.time_ms, df.S_m,  color=col, ls='--', lw=2,
            label=f'{grp} – gaze S')
    ax.fill_between(df.time_ms,
                    df.S_m-df.S_sem, df.S_m+df.S_sem,
                    color=col, alpha=.20)

ax.axhline(.5, ls='--', color='k')
ax.axvline(0,  ls='--', color='k')
ax.set(xlim=(0, rt_glob_mean+100), ylim=(0,1),
       xlabel='Time (ms)', ylabel='Fixation ratio',
       title='ES phase – gaze to E vs S, split by accuracy')
ax.legend(frameon=False, fontsize=7, ncol=1,
          bbox_to_anchor=(1.02,1), loc='upper left')
fig.tight_layout()
fig.savefig(os.path.join(fig_dir, 'ES_EvS_byAccuracy.png'),
            dpi=300, bbox_inches='tight')
plt.show()


# additonal plots per OV block -------------------------------------------

OV_level = {1: "low",
            2: "medium",
            3: "high"}  

block_cols = ['pink', 'violet', 'indigo']

for OV_val, col in zip(range(1, 4), block_cols):
    label = OV_level[OV_val]
    blk_raw = raw[raw['OV_2'] == OV_val] 
    if blk_raw.empty:
        print(f'Block {OV_val} is empty – skipped.')
        continue

    # choice-phase bins  ------------------------------------
    tbl = blk_raw.apply(trial_bins, axis=1).dropna()
    tbl = pd.concat(tbl.tolist(), ignore_index=True)

    rts = (blk_raw[['sub_id','TrialID','rtime']]
           .rename(columns={'sub_id':'sub','TrialID':'trial','rtime':'rt'}))
    rts['rt_bin'] = (rts.rt*1000 // BIN_MS).astype(int)
    tbl = tbl.merge(rts, on=['sub','trial'])
    tbl = tbl[tbl.bin <= tbl.rt_bin]
    tbl['time_ms'] = tbl['bin'] * BIN_MS

    sub_traj = (tbl.groupby(['sub','time_ms']).hi_ratio
                    .mean().reset_index(name='hi_ratio'))
    grand = (sub_traj.groupby('time_ms')
                     .agg(mean=('hi_ratio','mean'),
                          sem =('hi_ratio', lambda x: x.std()/np.sqrt(len(x))))
                     .reset_index())
    grand['lo_mean'] = 1 - grand['mean']
    grand['lo_sem']  = grand['sem']

    rt_mean = blk_raw.rtime.mean()*1000
    rt_sem  = blk_raw.rtime.std(ddof=1)/np.sqrt(len(blk_raw))*1000

    # for E vs S ######################################################################
    ES_tbl = blk_raw.apply(trial_bins_ES, axis=1).dropna()
    ES_tbl = pd.concat(ES_tbl.tolist(), ignore_index=True)

    rts_ES = (raw[['sub_id','TrialID','rtime']]
              .rename(columns={'sub_id':'sub','TrialID':'trial','rtime':'rt'}))
    rts_ES['rtES_bin'] = (rts_ES.rt*1000 // BIN_MS).astype(int)

    ES_tbl = ES_tbl.merge(rts_ES, on=['sub','trial'])
    ES_tbl = ES_tbl[ES_tbl.bin <= ES_tbl.rtES_bin]          

    # time in ms relative to onset
    ES_tbl['timeES_ms'] = ES_tbl['bin'] * BIN_MS

    sub_trajES = (ES_tbl.groupby(['sub','timeES_ms']).ES_ratio.mean().reset_index(name='ES_ratio'))

    grandES = (sub_trajES.groupby('timeES_ms')
                     .agg(mean=('ES_ratio','mean'),
                          sem =('ES_ratio', lambda x: x.std()/np.sqrt(len(x))))
                     .reset_index())
    grandES['S_mean'] = 1 - grandES['mean']
    grandES['S_sem']  = grandES['sem']

    rtES_mean = blk_raw.rtime.mean()*1000
    rtES_sem  = blk_raw.rtime.std(ddof=1)/np.sqrt(len(blk_raw))*1000

    # ----------------------  PLOTS  -----------------------------

    # choice figure -----------------------------------
    fig_c, ax_c = plt.subplots(figsize=(4.3,3))
    ax_c.axvspan(rt_mean-rt_sem, rt_mean+rt_sem,
                 color='lightgrey', alpha=.5, zorder=0)
    ax_c.axvline(rt_mean, color='grey', lw=1.2)

    ax_c.plot(grand.time_ms, grand['mean'],  color=col, lw=2)
    ax_c.fill_between(grand.time_ms,
                      grand['mean']-grand['sem'],
                      grand['mean']+grand['sem'],
                      color=col, alpha=.25)
    ax_c.plot(grand.time_ms, grand['lo_mean'], color=col,
              ls='--', lw=2)
    ax_c.fill_between(grand.time_ms,
                      grand['lo_mean']-grand['lo_sem'],
                      grand['lo_mean']+grand['lo_sem'],
                      color=col, alpha=.25)

    ax_c.axhline(.5, ls='--', color='k')
    ax_c.set(xlim=(0, rt_mean+100), ylim=(0,1),
             xlabel='Time (ms)', ylabel='Fixation ratio',
             title=f'OV-level {label}: Choice')

    # ES  E VS S figure -------------------
    fig_f, ax_f = plt.subplots(figsize=(4.3,3))
    ax_f.axvspan(rtES_mean-rtES_sem, rtES_mean+rtES_sem,
                 color='lightgrey', alpha=.5, zorder=0)
    ax_f.axvline(rtES_mean, color='grey', lw=1.2)    
    ax_f.plot(grandES.timeES_ms, grandES['mean'],  color=col, lw=2)
    
    ax_f.fill_between(grandES.timeES_ms,
                      grandES['mean']-grandES['sem'],
                      grandES['mean']+grandES['sem'],
                      color=col, alpha=.25)
    ax_f.plot(grandES.timeES_ms, grandES['S_mean'], color=col,
              ls='--', lw=2)
    ax_f.fill_between(grandES.timeES_ms,
                      grandES['S_mean']-grandES['S_sem'],
                      grandES['S_mean']+grandES['S_sem'],
                      color=col, alpha=.25)

    ax_f.axhline(.5, ls='--', color='k')
    ax_f.set(xlim=(0, rtES_mean+100), ylim=(0,1),
             xlabel='Time (ms)', ylabel='Fixation ratio',
             title=f'OV-level {label}: Choice')
    

    # legends
    solid_hi  = Line2D([], [], color=col, lw=2,      label='higher EV')
    dashed_lo = Line2D([], [], color=col, ls='--', lw=2, label='lower EV')

    solid_E = Line2D([], [], color=col, lw=2,      label='E')
    dash_S   = Line2D([], [], color=col, ls='--', lw=2, label='S')

    ax_c.legend(handles=[solid_hi, dashed_lo], frameon=False, fontsize=7)
    ax_f.legend(handles=[solid_E, dash_S], frameon=False, fontsize=7)

    fig_c.savefig(os.path.join(fig_dir, f'OV{label}_ES_EV.png'),dpi=300, bbox_inches='tight')
    fig_f.savefig(os.path.join(fig_dir, f'OV{label}_ES_ES.png'),dpi=300, bbox_inches='tight')




# -------------------------------------------------------------
all_curves = defaultdict(list)          # keys for choice,feed, outcf


for OV_val, col in zip(range(1, 4), block_cols):      
    label   = OV_level[OV_val]
    blk_raw = raw[raw['OV_2'] == OV_val]             
    if blk_raw.empty:
        continue
    # choice: hi vs lo EV -------------------------
    tbl   = blk_raw.apply(trial_bins, axis=1).dropna()
    tbl   = pd.concat(tbl.tolist(), ignore_index=True)
    rts   = blk_raw[['sub_id','TrialID','rtime']].rename(
              columns={'sub_id':'sub','TrialID':'trial','rtime':'rt'})
    rts['rt_bin'] = (rts.rt*1000 // BIN_MS).astype(int)
    tbl   = tbl.merge(rts, on=['sub','trial'])
    tbl   = tbl[tbl.bin <= tbl.rt_bin]
    tbl['time_ms'] = tbl['bin']*BIN_MS
    sub   = tbl.groupby(['sub','time_ms']).hi_ratio.mean().reset_index()
    g     = sub.groupby('time_ms').hi_ratio.agg(['mean','sem']).reset_index()
    g['lo_mean'] = 1-g['mean'];  g['lo_sem'] = g['sem']
    all_curves['choice'].append((label, col, g, 'higher EV', 'lower EV'))

    # feedback: hi vs lo EV -----------------------
    
    tbl_ES   = blk_raw.apply(trial_bins_ES, axis=1).dropna()
    tbl_ES   = pd.concat(tbl_ES.tolist(), ignore_index=True)

    rts_ES   = blk_raw[['sub_id','TrialID','rtime']].rename(
              columns={'sub_id':'sub',
                       'TrialID':'trial',
                       'rtime':'rt'})
    rts_ES['rtES_bin'] = (rts_ES.rt * 1000 // BIN_MS).astype(int)

    tbl_ES = tbl_ES.merge(rts_ES, on=['sub','trial'])
    tbl_ES = tbl_ES[tbl_ES.bin <= tbl_ES.rtES_bin]   
    tbl_ES['timeES_ms'] = tbl_ES['bin'] * BIN_MS
    sub_ES   = tbl_ES.groupby(['sub','timeES_ms']).ES_ratio.mean().reset_index()
    g_ES     = sub_ES.groupby('timeES_ms').ES_ratio.agg(['mean','sem']).reset_index()
    g_ES['S_mean'] = 1-g_ES['mean'];  g_ES['S_sem'] = g_ES['sem']
    g_ES.rename(columns={'timeES_ms': 'time_ms'}, inplace=True)  
    all_curves['ES'].append((label, col, g_ES, 'E', 'S'))   

# helper to draw one overlay figure
# ---------------------------------------------------------------
def overlay_plot(key, title, ylab, fname):
    
    #key        : choice, ES
    #title      : figure title
    #ylab       : y-axis label
    #fname      : file name (in fig_dir)
    
    fig, ax = plt.subplots(figsize=(6.2, 3.8))
    
    ax.axvspan(rt_glob_mean-rt_glob_sem,
               rt_glob_mean+rt_glob_sem,
               color='lightgrey', alpha=.35, zorder=10)
    ax.axvline(rt_glob_mean, color='grey', lw=1.2, zorder=11)

    for label, col, df, solid_lbl, dash_lbl in all_curves[key]:
        # solid line + SEM band
        ax.plot(df.time_ms, df['mean'],
                color=col, lw=2, label=f'{label} – {solid_lbl}')
        ax.fill_between(df.time_ms,
                        df['mean']-df['sem'],
                        df['mean']+df['sem'],
                        color=col, alpha=.20)

        # dashed line + SEM band
        second   = 'lo_mean' if key == 'choice' else 'S_mean'
        second_s = 'lo_sem'  if key == 'choice' else 'S_sem'
        ax.plot(df.time_ms, df[second],
                color=col, ls='--', lw=2, label=f'{label} – {dash_lbl}')
        ax.fill_between(df.time_ms,
                        df[second]-df[second_s],
                        df[second]+df[second_s],
                        color=col, alpha=.20)

    ax.axhline(.5, ls='--', color='k')
    ax.axvline(0 , ls='--', color='k')
    ax.set(xlim=(0, 1300), ylim=(0,1),
           xlabel='Time (ms)', ylabel=ylab, title=title)

    ax.legend(frameon=False, fontsize=7, ncol=1,
              bbox_to_anchor=(1.02, 1), loc="upper left", borderaxespad=0.)
    fig.tight_layout()
    fig.savefig(os.path.join(fig_dir, fname),
                dpi=300, bbox_inches='tight')
    plt.show()

#  three overlay figures
# ---------------------------------------------------------------
overlay_plot('choice',
             'Choice – gaze on higher vs lower EV (OV levels)',
             'Fixation ratio',
             'ALL_OV_levels_ES_EV.png')

overlay_plot('ES',
             'Choice – gaze on E vs S (OV levels)',
             'Fixation ratio',
             'ALL_OV_levels_ES_ES.png')


#  FIXATION-BY-FIXATION (response-aligned, ES phase only)


EXCLUDE_SUBS = {1, 4, 5, 6, 14, 99}
if "sub_id" not in raw.columns:
    raise KeyError("Expected column sub_id in raw.")
raw = raw[~raw["sub_id"].isin(EXCLUDE_SUBS)].copy()
print(f"Excluded subjects: {sorted(EXCLUDE_SUBS)}")
print(f"Included subjects: {sorted(raw['sub_id'].unique().tolist())}")

# helper
def gather_fixseq(row, back=4):
    # list with first element (index 0) beign the fixation at the moment of response
    # all fixtions are fixations preceding it row.fixations is a list of dicts or a string
    # that can be parsed. Each dict contains a start_time in ms relative to first fix, roi which
    # is either E=1 or S=2 & row.time = response time in sec
    
    fx_raw = row.Fixations
    if pd.isna(fx_raw):
        return []

    fixes = fx_raw if isinstance(fx_raw, list) else ast.literal_eval(fx_raw)
    if not fixes:
        return []

    t0    = fixes[0]['start_time']
    rt_ms = row.rtime * 1000 

    usable = [f for f in fixes if (f['start_time'] - t0) <= rt_ms]
    if not usable:
        return []

    lab = {1: "E", 2: "S"}
    seq = [lab[f['roi']] for f in usable if f.get('roi') in (1, 2)]
    seq = seq[::-1]
    return seq[:back + 1]

#  build sequences & drop trials with <1 fixation (but this can actually be changed)

MAX_BACK = 4
raw_ES            = raw.copy()
raw_ES['fix_seq'] = raw_ES.apply(lambda r: gather_fixseq(r, back=MAX_BACK), axis=1)
raw_ES = raw_ES[raw_ES.fix_seq.apply(len) > -1]
if raw_ES.empty:
    raise RuntimeError("No ES trial has usable fixations after exclusions")

# long table: one row = one fixation

records = []
for _, r in raw_ES.iterrows():
    for idx, lab in enumerate(r.fix_seq):       
        if idx > MAX_BACK:                   
            break
        records.append(dict(sub     = r.sub_id,
                            fix_idx = idx,      # 0, 1, 2, 3, 4
                            is_E    = int(lab == "E"),
                            is_S    = int(lab == "S")))
fix_df = pd.DataFrame.from_records(records)

# participant means that is grand mean & SEM
# ------------------------------------------------------------------------
sub_fix = (fix_df.groupby(["sub", "fix_idx"], as_index=False)
                   .agg(E_ratio=("is_E", "mean"),
                        S_ratio=("is_S", "mean")))

def _sem(x):
    x = pd.Series(x)
    return x.std(ddof=1) / np.sqrt(len(x)) if len(x) > 1 else np.nan

grand_fix = (sub_fix.groupby("fix_idx", as_index=False)
                      .agg(E_m=("E_ratio", "mean"),
                           E_sem=("E_ratio", _sem),
                           S_m=("S_ratio", "mean"),
                           S_sem=("S_ratio", _sem)))

# keep only indices up to MAX_BACK (0..4), then convert to 0, -1, -2, ...
grand_fix = grand_fix[grand_fix["fix_idx"] <= MAX_BACK].copy()
grand_fix["x"] = -grand_fix.fix_idx
grand_fix.loc[grand_fix.fix_idx == 0, "x"] = 0

# plot
# ------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(4.8, 3.0))

ax.plot(grand_fix.x, grand_fix.E_m, color="deepskyblue", lw=2, label="E")
ax.plot(grand_fix.x, grand_fix.S_m, color="darkorchid", ls="--", lw=2, label="S")

ax.fill_between(grand_fix.x,
                grand_fix.E_m - grand_fix.E_sem,
                grand_fix.E_m + grand_fix.E_sem,
                color="deepskyblue", alpha=.25)
ax.fill_between(grand_fix.x,
                grand_fix.S_m - grand_fix.S_sem,
                grand_fix.S_m + grand_fix.S_sem,
                color="darkorchid", alpha=.25)

ax.axhline(.5, ls="--", color="k")
ax.axvline(0,  color="grey", lw=1.2)

ax.set(xlim=(-MAX_BACK, .1), ylim=(0, 1),
       xlabel="Fixation index (0 = at response, negative = before)",
       ylabel="Probability of fixating option",
       title=f"ES phase – fixation at response and previous {MAX_BACK}")

ax.set_xticks([0] + [-i for i in range(1, MAX_BACK + 1)])
ax.legend(frameon=False, fontsize=8, title="Gaze on")
fig.tight_layout()

fig_fixseq = os.path.join(fig_dir, f"ES_fixseq_clickPlus{MAX_BACK}_noSig_OV.png")
fig.savefig(fig_fixseq, dpi=300, bbox_inches="tight")
plt.show()
print("Saved:", fig_fixseq)

# paired t-tests is E_ratio not S_ratio across subjects per fixation idx?
# ------------------------------------------------------------------------
ttest_results = []
for i in range(MAX_BACK + 1):
    data_i = sub_fix[sub_fix.fix_idx == i]
    if data_i.empty:
        ttest_results.append((np.nan, np.nan))
        continue
    t_stat, p_val = ttest_rel(data_i.E_ratio, data_i.S_ratio)
    ttest_results.append((t_stat, p_val))

# Add t/p-values to grand_fix
t_vals, p_vals = zip(*ttest_results)
grand_fix["t_value"] = t_vals
grand_fix["p_value_uncorrected"] = p_vals

# FDR correction
_, p_vals_corr, _, _ = multipletests(p_vals, method="fdr_bh")
grand_fix["p_value_fdr"] = p_vals_corr

# export table (mean & SEM + t/p values)
# ------------------------------------------------------------------------
n_subs = sub_fix["sub"].nunique()
grand_table = (grand_fix
               .loc[:, ["x", "E_m", "E_sem", "S_m", "S_sem",
                        "t_value", "p_value_uncorrected", "p_value_fdr"]]
               .assign(n_subs=n_subs)
               .sort_values("x"))

grand_csv = os.path.join(fig_dir, f"ES_fixseq_grand_mean_clickPlus{MAX_BACK}_withSig_OV.csv")
grand_table.to_csv(grand_csv, index=False)
print(f"Grand-mean fixation table with stats saved in {grand_csv}")