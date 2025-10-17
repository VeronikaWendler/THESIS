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



# big file that has everything
df_raw = pd.read_csv("D:/Aberdeen_Uni_June24/cap/THESIS/Garcia_Analysis/data/data_sets/GarciaParticipants_Eye_Response_Feed_Allfix_addm_OV_Abs_CCT.csv")

# figures dir
fig_dir = "D:/Aberdeen_Uni_June24/cap/THESIS/Garcia_Analysis/Figures/gaze_trajectory"


BIN_MS    = 25                       # bin width in ms
ROI   = {1: 'left', 2: 'right'}      # ROI

raw = (df_raw
         .query("phase == 'ES' and p1 != 50")      
         .copy())

# pie_split = 1  (p2  > 50 %)   |  0  (p2  < 50 %)
raw['Im_split'] = (raw['p1'] > 50).astype(int)

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
    """
    Return per-25 ms bins:
        • E_cnt, S_cnt
        • identify trial group: 'E-chosen' (cho==1) or 'S-chosen' (cho==2)
    """
    try:
        fixes = ast.literal_eval(row.Fixations) or []
    except Exception:
        return None
    if not fixes:
        return None

    E_side, S_side = 'left', 'right'
    t0, bins = fixes[0]['start_time'], {}

    for f in fixes:
        side = ROI.get(f['roi'], 'none')            # 'left' | 'right' | 'none'
        idx  = 0 if side == E_side else 1           # 0 = E, 1 = S
        for b in range(f['start_time']//BIN_MS - t0//BIN_MS,
                       (f['end_time']-1)//BIN_MS - t0//BIN_MS + 1):
            if b not in bins:
                bins[b] = [0, 0]                    # [E_cnt, S_cnt]
            bins[b][idx] += 1

    out = pd.DataFrame.from_dict(bins, orient='index',
                                 columns=['E_cnt', 'S_cnt'])
    out['sub']      = row.sub_id
    out['trial']    = row.TrialID
    out['bin']      = out.index
    out['cho_grp']  = 'E-chosen' if row.cho == 1 else 'S-chosen'
    out['OV_lvl']  = row.OV_2          
    out['VD_lvl']  = row.VD_2          
    out['ratio_E']  = out.E_cnt / (out.E_cnt + out.S_cnt)
    out['ratio_S']  = 1 - out['ratio_E']
    return out[['sub','trial','bin','OV_lvl','cho_grp','VD_lvl','ratio_E','ratio_S']]




def trial_bins_EvS_byAcc(row):
    """
    For every 25-ms bin during ES trials count E- and S-fixations.
    Attach a trial label:
        • 'Correct'    (row.corr == 1)
        • 'Incorrect'  (row.corr == 0)
    """
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
    out['OV_lvl']  = row.OV_2
    out['VD_lvl']  = row.VD_2
    out['ratio_E'] = out.E_cnt / (out.E_cnt + out.S_cnt)
    out['ratio_S'] = 1 - out['ratio_E']
    return out[['sub','trial','bin','acc_grp','OV_lvl','VD_lvl','ratio_E','ratio_S']]



def trial_bins_EvS_byIM(row):
    """
    For every 25-ms bin count fixations on E vs S
    and tag the trial as pie>0.5 (1) or pie<0.5 (0).
    """

    try:
        fixes = ast.literal_eval(row.Fixations) or []
    except Exception:
        return None
    if not fixes:
        return None

    E_side, S_side = 'left', 'right'
    t0, bins = fixes[0]['start_time'], {}

    for f in fixes:
        idx = 0 if ROI.get(f['roi'],'none') == E_side else 1   # 0 = E, 1 = S
        for b in range(f['start_time']//BIN_MS - t0//BIN_MS,
                       (f['end_time']-1)//BIN_MS - t0//BIN_MS + 1):
            bins.setdefault(b, [0,0])[idx] += 1

    out = pd.DataFrame.from_dict(bins, orient='index',
                                 columns=['E_cnt','S_cnt'])
    out['sub']       = row.sub_id
    out['trial']     = row.TrialID
    out['bin']       = out.index
    out['IM_grp']   = 'Image>0.5' if row.Im_split == 1 else 'Image<0.5'
    out['OV_lvl']    = row.OV_2
    out['VD_lvl']  = row.VD_2
    out['ratio_E']   = out.E_cnt / (out.E_cnt + out.S_cnt)
    out['ratio_S']   = 1 - out['ratio_E']
    return out[['sub','trial','bin','IM_grp','OV_lvl','VD_lvl','ratio_E','ratio_S']]



# bin table per trial -------------------------------------------------------------------------------------
tbl = raw.apply(trial_bins, axis=1).dropna()
tbl = pd.concat(tbl.tolist(), ignore_index=True)

# mask bins after each individual RT 
rts = (raw[['sub_id','TrialID','rtime']]
       .rename(columns={'sub_id':'sub','TrialID':'trial','rtime':'rt'}))
rts['rt_bin'] = (rts.rt*1000 // BIN_MS).astype(int)

tbl = tbl.merge(rts, on=['sub','trial'])
tbl = tbl[tbl.bin <= tbl.rt_bin]           # keep up to the click

# time in ms relative to onset
tbl['time_ms'] = tbl['bin'] * BIN_MS

# subject-average trajectory 
sub_traj = (tbl.groupby(['sub','time_ms']).hi_ratio.mean().reset_index(name='hi_ratio'))

# grand mean & SEM across participants
grand = (sub_traj.groupby('time_ms')
                   .agg(mean=('hi_ratio','mean'),
                        sem =('hi_ratio', lambda x: x.std()/np.sqrt(len(x))))
                   .reset_index())
grand['lo_mean'] = 1 - grand['mean']   # complementary curve
grand['lo_sem' ] = grand['sem']        # same SEM

# global mean RT & its SEM 
rt_mean = raw.rtime.mean()*1000            # ms
rt_sem  = raw.rtime.std(ddof=1)/np.sqrt(len(raw))*1000




# ─── build table for ES phase ─────────────────────────────────────
ES_tbl = raw.apply(trial_bins_ES, axis=1).dropna()
ES_tbl = pd.concat(ES_tbl.tolist(), ignore_index=True)

rts_ES = (raw[['sub_id','TrialID','rtime']]
       .rename(columns={'sub_id':'sub','TrialID':'trial','rtime':'rt'}))
rts_ES['rtES_bin'] = (rts_ES.rt*1000 // BIN_MS).astype(int)

ES_tbl = ES_tbl.merge(rts_ES, on=['sub','trial'])
ES_tbl = ES_tbl[ES_tbl.bin <= ES_tbl.rtES_bin]           # keep up to the click

# time in ms relative to onset
ES_tbl['timeES_ms'] = ES_tbl['bin'] * BIN_MS

# subject-average trajectory 
sub_trajES = (ES_tbl.groupby(['sub','timeES_ms']).ES_ratio.mean().reset_index(name='ES_ratio'))

# grand mean & SEM across participants
grandES = (sub_trajES.groupby('timeES_ms')
                   .agg(mean=('ES_ratio','mean'),
                        sem =('ES_ratio', lambda x: x.std()/np.sqrt(len(x))))
                   .reset_index())
grandES['S_mean'] = 1 - grandES['mean']   # complementary curve
grandES['S_sem' ] = grandES['sem']        # same SEM

# global mean RT & its SEM 
rtES_mean = raw.rtime.mean()*1000            # ms
rtES_sem  = raw.rtime.std(ddof=1)/np.sqrt(len(raw))*1000



#----------------------------------------------------------------------------
# chocie E vs S

t_cap = int(rt_glob_mean + 100)      # ≈ 1 × mean RT + 100 ms
choice_tbl = raw.apply(trial_bins_EvS_byChoice, axis=1).dropna()
choice_tbl = pd.concat(choice_tbl.tolist(), ignore_index=True)
choice_tbl['time_ms'] = choice_tbl['bin'] * BIN_MS
choice_tbl = choice_tbl[choice_tbl.time_ms <= rt_glob_mean + 100]  # trim

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



# ----------  accuracy-based E vs S  ---------- #
acc_tbl = raw.apply(trial_bins_EvS_byAcc, axis=1).dropna()
acc_tbl = pd.concat(acc_tbl.tolist(), ignore_index=True)
acc_tbl['time_ms'] = acc_tbl['bin'] * BIN_MS
acc_tbl = acc_tbl[acc_tbl.time_ms <= rt_glob_mean + 100]   # trim to mean RT + 100 ms

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



# ----- build table & trim to RT+100 ms ------------------------------------------
IM_tbl = raw.apply(trial_bins_EvS_byIM, axis=1).dropna()
IM_tbl = pd.concat(IM_tbl.tolist(), ignore_index=True)
IM_tbl['time_ms'] = IM_tbl['bin'] * BIN_MS
IM_tbl = IM_tbl[IM_tbl.time_ms <= rt_glob_mean + 100]

# subject-level
sub_IM = (IM_tbl
           .groupby(['sub','IM_grp','time_ms'])
           .agg(E_mean=('ratio_E','mean'),
                S_mean=('ratio_S','mean'))
           .reset_index())

# grand mean + SEM
grp_IM = (sub_IM
           .groupby(['IM_grp','time_ms'])
           .agg(E_m=('E_mean','mean'),
                E_sem=('E_mean', lambda x: x.std(ddof=1)/np.sqrt(len(x))),
                S_m=('S_mean','mean'),
                S_sem=('S_mean', lambda x: x.std(ddof=1)/np.sqrt(len(x))))
           .reset_index())


# -----------------------------------------------------------------------
# helper: small, tidy legend kwargs
legend_kw = dict(loc='upper right', frameon=False,
                 fontsize=8,  title_fontsize=9,
                 borderpad=0.2, labelspacing=0.25, handlelength=2)

line_handles = [
    Line2D([], [], color='darkorchid',  lw=2, label='higher'),
    Line2D([], [], color='darkorchid', ls='--', lw=2, label='lower')
]

line_handles2 = [
    Line2D([], [], color='royalblue',  lw=2, label='E'),
    Line2D([], [], color='darkorchid', ls='--', lw=2, label='S')
]

# ────────────────────────────────────────────────────────────────────────
# E or S chosen EV
# ────────────────────────────────────────────────────────────────────────
fig_choice, ax_c = plt.subplots(figsize=(4.5, 3))

# RT band  (mean ± SEM)
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

# ────────────────────────────────────────────────────────────────────────
# FEEDBACK  (outcome-locked) – separate figure
# ────────────────────────────────────────────────────────────────────────

fig_choiceES, ax_ES = plt.subplots(figsize=(4.5, 3))

# RT band  (mean ± SEM)
ax_ES.axvspan(rtES_mean-rtES_sem, rtES_mean+rtES_sem,
             color='lightgrey', alpha=.5, zorder=0)
# RT mean line on top of the band
ax_ES.axvline(rtES_mean, color='grey', lw=1.2, zorder=1)

# higher-EV
ax_ES.plot(grandES.timeES_ms, grandES['mean'],
          color='royalblue', lw=2)
ax_ES.fill_between(grandES.timeES_ms,
                  grandES['mean']-grandES['sem'],
                  grandES['mean']+grandES['sem'],
                  color='royalblue', alpha=.25, label='_nolegend_')

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

ax_ES.legend(handles=line_handles2, title='Gaze on', **legend_kw)
fig_choiceES.tight_layout()
plot_path2= os.path.join(fig_dir, "gaze_choice_ES_ES.png")
fig_choiceES.savefig(plot_path2, dpi=300, bbox_inches='tight')
fig_choiceES.show()

#-----------------------------------------------------------------------------------------------------------------
# E vs S by chosen item

fig, ax = plt.subplots(figsize=(4.8,3.3))

# RT cue
ax.axvspan(rt_glob_mean-rt_glob_sem, rt_glob_mean+rt_glob_sem,
           color='lightgrey', alpha=.35, zorder=0)
ax.axvline(rt_glob_mean, color='grey', lw=1.2, zorder=1)

colour_grp = {'E-chosen':'royalblue', 'S-chosen':'darkorange'}

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
# E vs S by accuracy 

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



# plot pie>05, vs pie<0.5
fig, ax = plt.subplots(figsize=(4.8,3.3))

ax.axvspan(rt_glob_mean-rt_glob_sem, rt_glob_mean+rt_glob_sem,
           color='lightgrey', alpha=.35, zorder=0)
ax.axvline(rt_glob_mean, color='grey', lw=1.2, zorder=1)

colours = {'Image>0.5':'mediumorchid', 'Image<0.5':'darkturquoise'}

for grp, df in grp_IM.groupby('IM_grp'):
    c = colours[grp]

    ax.plot(df.time_ms, df.E_m,  color=c, lw=2,
            label=f'{grp} – gaze E')
    ax.fill_between(df.time_ms, df.E_m-df.E_sem, df.E_m+df.E_sem,
                    color=c, alpha=.20)

    ax.plot(df.time_ms, df.S_m,  color=c, ls='--', lw=2,
            label=f'{grp} – gaze S')
    ax.fill_between(df.time_ms, df.S_m-df.S_sem, df.S_m+df.S_sem,
                    color=c, alpha=.20)

ax.axhline(.5, ls='--', color='k')
ax.axvline(0,  ls='--', color='k')
ax.set(xlim=(0, rt_glob_mean+100), ylim=(0,1),
       xlabel='Time (ms)', ylabel='Fixation ratio',
       title='ES phase – gaze E vs S, split by Image-probability')
ax.legend(frameon=False, fontsize=7, ncol=1,
          bbox_to_anchor=(1.02,1), loc='upper left')
fig.tight_layout()
fig.savefig(os.path.join(fig_dir, 'ES_EvS_byImage.png'),
            dpi=300, bbox_inches='tight')
plt.show()



# plots per LE block

# ──────────────────────────────────────────────────────────────────────────
#  EXTRA: plots *per learning block*  (cond = 0‒3)
# ──────────────────────────────────────────────────────────────────────────

# ---------------------------------------------------------------
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

    # ===========================================================
    # ---------- RE-RUN the pipeline on *blk_raw*  ---------------
    # ===========================================================
    # 1)  choice-phase bins  ------------------------------------
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
    ES_tbl = ES_tbl[ES_tbl.bin <= ES_tbl.rtES_bin]           # keep up to the click

    # time in ms relative to onset
    ES_tbl['timeES_ms'] = ES_tbl['bin'] * BIN_MS

    # subject-average trajectory 
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
    # ===========================================================

    # ----- figure 1: choice -----------------------------------
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

    # ----- figure 2: ES  E VS S  -------------------
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
    
    # tidy legends (same style for every plot)
    # ── tidy legends (same style for every plot) ─────────────────────────────
    solid_hi  = Line2D([], [], color=col, lw=2,      label='higher EV')
    dashed_lo = Line2D([], [], color=col, ls='--', lw=2, label='lower EV')

    solid_E = Line2D([], [], color=col, lw=2,      label='E')
    dash_S   = Line2D([], [], color=col, ls='--', lw=2, label='S')

    ax_c.legend(handles=[solid_hi, dashed_lo], frameon=False, fontsize=7)
    ax_f.legend(handles=[solid_E, dash_S], frameon=False, fontsize=7)

    fig_c.savefig(os.path.join(fig_dir, f'OV{label}_ES_EV.png'),dpi=300, bbox_inches='tight')
    fig_f.savefig(os.path.join(fig_dir, f'OV{label}_ES_ES.png'),dpi=300, bbox_inches='tight')




# 1) -------------------------------------------------------------
all_curves = defaultdict(list)          # keys → 'choice' | 'feed' | 'outcf'


for OV_val, col in zip(range(1, 4), block_cols):      # 1, 2, 3  (= low-medium-high)
    label   = OV_level[OV_val]
    blk_raw = raw[raw['OV_2'] == OV_val]              # ← FIX 1
    if blk_raw.empty:
        continue
    # ---------- a) choice: hi vs lo EV -------------------------
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

    # ---------- b) feedback: hi vs lo EV -----------------------
    
    tbl_ES   = blk_raw.apply(trial_bins_ES, axis=1).dropna()
    tbl_ES   = pd.concat(tbl_ES.tolist(), ignore_index=True)

    rts_ES   = blk_raw[['sub_id','TrialID','rtime']].rename(
              columns={'sub_id':'sub',
                       'TrialID':'trial',
                       'rtime':'rt'})
    rts_ES['rtES_bin'] = (rts_ES.rt * 1000 // BIN_MS).astype(int)

    tbl_ES = tbl_ES.merge(rts_ES, on=['sub','trial'])
    tbl_ES = tbl_ES[tbl_ES.bin <= tbl_ES.rtES_bin]   # ← now this works
    tbl_ES['timeES_ms'] = tbl_ES['bin'] * BIN_MS
    sub_ES   = tbl_ES.groupby(['sub','timeES_ms']).ES_ratio.mean().reset_index()
    g_ES     = sub_ES.groupby('timeES_ms').ES_ratio.agg(['mean','sem']).reset_index()
    g_ES['S_mean'] = 1-g_ES['mean'];  g_ES['S_sem'] = g_ES['sem']
    g_ES.rename(columns={'timeES_ms': 'time_ms'}, inplace=True)   # ← unify x-axis name
    all_curves['ES'].append((label, col, g_ES, 'E', 'S'))   


# ---------------------------------------------------------------
# 2)  helper to draw one overlay figure
# ---------------------------------------------------------------
def overlay_plot(key, title, ylab, fname):
    """
    key        : 'choice', 'ES'
    title      : figure title
    ylab       : y-axis label
    fname      : file name (saved in fig_dir)
    """
    fig, ax = plt.subplots(figsize=(6.2, 3.8))
    
    ax.axvspan(rt_glob_mean-rt_glob_sem,
               rt_glob_mean+rt_glob_sem,
               color='lightgrey', alpha=.35, zorder=10)
    ax.axvline(rt_glob_mean, color='grey', lw=1.2, zorder=11)

    for label, col, df, solid_lbl, dash_lbl in all_curves[key]:
        # solid line (+ SEM band)
        ax.plot(df.time_ms, df['mean'],
                color=col, lw=2, label=f'{label} – {solid_lbl}')
        ax.fill_between(df.time_ms,
                        df['mean']-df['sem'],
                        df['mean']+df['sem'],
                        color=col, alpha=.20)

        # dashed line (+ SEM band)
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

    # put legend outside (to the right)
    ax.legend(frameon=False, fontsize=7, ncol=1,
              bbox_to_anchor=(1.02, 1), loc="upper left", borderaxespad=0.)
    fig.tight_layout()
    fig.savefig(os.path.join(fig_dir, fname),
                dpi=300, bbox_inches='tight')
    plt.show()

# ---------------------------------------------------------------
# 3)  make the three overlay figures
# ---------------------------------------------------------------
overlay_plot('choice',
             'Choice – gaze on higher vs lower EV (OV levels)',
             'Fixation ratio',
             'ALL_OV_levels_ES_EV.png')

overlay_plot('ES',
             'Choice – gaze on E vs S (OV levels)',
             'Fixation ratio',
             'ALL_OV_levels_ES_ES.png')



########################################################################################################################################

# for chosen item #######################################################################
    
choice_tbl  = raw.apply(trial_bins_EvS_byChoice, axis=1).dropna()
choice_tbl  = pd.concat(choice_tbl.tolist(), ignore_index=True)
choice_tbl['time_ms'] = choice_tbl['bin'] * BIN_MS
choice_tbl  = choice_tbl[choice_tbl.time_ms <= rt_glob_mean + 100]

sub_choice = (choice_tbl
              .groupby(['sub','OV_lvl','cho_grp','time_ms'])
              .agg(E_mean=('ratio_E','mean'),
                   S_mean=('ratio_S','mean'))
              .reset_index())

grp_choice_OV = (sub_choice
                 .groupby(['OV_lvl','cho_grp','time_ms'])
                 .agg(E_m=('E_mean','mean'),
                      E_sem=('E_mean', lambda x: x.std(ddof=1)/np.sqrt(len(x))),
                      S_m=('S_mean','mean'),
                      S_sem=('S_mean', lambda x: x.std(ddof=1)/np.sqrt(len(x))))
                 .reset_index()) 
    

# for accuracy #######################################################################
    
    
acc_tbl   = raw.apply(trial_bins_EvS_byAcc, axis=1).dropna()
acc_tbl   = pd.concat(acc_tbl .tolist(), ignore_index=True)
acc_tbl ['time_ms'] = acc_tbl ['bin'] * BIN_MS
acc_tbl   = acc_tbl [acc_tbl .time_ms <= rt_glob_mean + 100]

sub_acc = (acc_tbl 
           .groupby(['sub','OV_lvl','acc_grp','time_ms'])
           .agg(E_mean=('ratio_E','mean'),
                S_mean=('ratio_S','mean'))
           .reset_index())

grp_acc_OV = (sub_acc
              .groupby(['OV_lvl','acc_grp','time_ms'])
              .agg(E_m=('E_mean','mean'),
                   E_sem=('E_mean', lambda x: x.std(ddof=1)/np.sqrt(len(x))),
                   S_m=('S_mean','mean'),
                   S_sem=('S_mean', lambda x: x.std(ddof=1)/np.sqrt(len(x))))
              .reset_index()) 
    
    
# for pie>0.5 pie<0.5 #######################################################################
    
IM_tbl   = raw.apply(trial_bins_EvS_byIM, axis=1).dropna()
IM_tbl   = pd.concat(IM_tbl .tolist(), ignore_index=True)
IM_tbl['time_ms'] = IM_tbl['bin'] * BIN_MS
IM_tbl   = IM_tbl[IM_tbl.time_ms <= rt_glob_mean + 100]

sub_IM = (IM_tbl 
           .groupby(['sub','OV_lvl','IM_grp','time_ms'])
           .agg(E_mean=('ratio_E','mean'),
                S_mean=('ratio_S','mean'))
           .reset_index())
grp_IM_OV = (sub_IM
              .groupby(['OV_lvl','IM_grp','time_ms'])
              .agg(E_m=('E_mean','mean'),
                   E_sem=('E_mean', lambda x: x.std(ddof=1)/np.sqrt(len(x))),
                   S_m=('S_mean','mean'),
                   S_sem=('S_mean', lambda x: x.std(ddof=1)/np.sqrt(len(x))))
              .reset_index()) 


# ===========================================







OV_cols = {1: 'pink', 2: 'violet', 3: 'indigo'}   #⇐ your colour map

def overlay_OV(df, split_col, title, fname,
               solid_tag, dashed_tag, solid_lab, dashed_lab):
    """
    • df         : grp_*_OV table (has OV_lvl, time_ms, E_m, S_m …)
    • split_col  : 'cho_grp' | 'acc_grp' | 'pie_grp'
    • solid_tag  : value in split_col that will be plotted solid (E-gaze)
    • dashed_tag : value in split_col that will be dashed  (S-gaze)
    """
    fig, ax = plt.subplots(figsize=(6.0,3.6))
    ax.axvspan(rt_glob_mean-rt_glob_sem, rt_glob_mean+rt_glob_sem,
               color='lightgrey', alpha=.35)
    ax.axvline(rt_glob_mean, color='grey', lw=1.2)

    for ov in (1,2,3):                                   # low-med-high
        col = OV_cols[ov]
        d   = df[df.OV_lvl == ov]

        # ── solid (E) ──────────────────────────────────────────────
        s = d[d[split_col] == solid_tag]
        ax.plot(s.time_ms, s.E_m,  color=col, lw=2,
                label=f'OV-{OV_level[ov]} – {solid_lab}')
        ax.fill_between(s.time_ms, s.E_m-s.E_sem, s.E_m+s.E_sem,
                        color=col, alpha=.20)

        # ── dashed (S) ─────────────────────────────────────────────
        dsh = d[d[split_col] == dashed_tag]
        ax.plot(dsh.time_ms, dsh.S_m, color=col, ls='--', lw=2,
                label=f'OV-{OV_level[ov]} – {dashed_lab}')
        ax.fill_between(dsh.time_ms, dsh.S_m-dsh.S_sem, dsh.S_m+dsh.S_sem,
                        color=col, alpha=.20)

    ax.axhline(.5, ls='--', color='k');  ax.axvline(0, ls='--', color='k')
    ax.set(xlim=(0, rt_glob_mean+100), ylim=(0,1),
           xlabel='Time (ms)', ylabel='Fixation ratio', title=title)
    ax.legend(frameon=False, fontsize=7, ncol=1,
              bbox_to_anchor=(1.02,1), loc='upper left')
    fig.tight_layout()
    fig.savefig(os.path.join(fig_dir, fname), dpi=300, bbox_inches='tight')
    plt.show()




# CHOSEN-item overlay
overlay_OV(grp_choice_OV, 'cho_grp',
           'Gaze E vs S – split by chosen item (OV overlay)',
           'OVoverlay_ChoiceSplit.png',
           solid_tag='E-chosen', dashed_tag='S-chosen',
           solid_lab='E-chosen | gaze E', dashed_lab='E-chosen | gaze S')

# ACCURACY overlay
overlay_OV(grp_acc_OV, 'acc_grp',
           'Gaze E vs S – split by accuracy (OV overlay)',
           'OVoverlay_AccuracySplit.png',
           solid_tag='Correct', dashed_tag='Incorrect',
           solid_lab='Correct | gaze E', dashed_lab='Correct | gaze S')

# PIE probability overlay
overlay_OV(grp_IM_OV, 'IM_grp',
           'Gaze E vs S – split by Image probability (OV overlay)',
           'OVoverlay_ImageSplit.png',
           solid_tag='Image>0.5', dashed_tag='Image<0.5',
           solid_lab='Image>0.5 | gaze E', dashed_lab='Image>0.5 | gaze S')
