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
df_raw = pd.read_csv("D:/Aberdeen_Uni_June24/cap/THESIS/OV_Analysis/data/data_sets/OVParticipants_Eye_Response_Feed_Allfix_addm_OV_Abs_CCT.csv")

# figures dir
fig_dir = "D:/Aberdeen_Uni_June24/cap/THESIS/OV_Analysis/Figures/gaze_trajectory"



BIN_MS    = 25                       # bin width in ms
ROI   = {1: 'left', 2: 'right'}      # ROI

raw = df_raw.query("phase == 'EE'").copy()   
rt_glob_mean = raw.rtime.mean() * 1000          
rt_glob_sem  = raw.rtime.std(ddof=1) / np.sqrt(len(raw)) * 1000

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



# -----------------------------------------------------------------------
# helper: small, tidy legend kwargs
legend_kw = dict(loc='upper right', frameon=False,
                 fontsize=8,  title_fontsize=9,
                 borderpad=0.2, labelspacing=0.25, handlelength=2)

line_handles = [
    Line2D([], [], color='royalblue',  lw=2, label='higher'),
    Line2D([], [], color='royalblue', ls='--', lw=2, label='lower')
]

# ────────────────────────────────────────────────────────────────────────
# CHOICE  (stimulus-locked) – its own figure
# ────────────────────────────────────────────────────────────────────────
fig_choice, ax_c = plt.subplots(figsize=(4.5, 3))

# RT band  (mean ± SEM)
ax_c.axvspan(rt_mean-rt_sem, rt_mean+rt_sem,
             color='lightgrey', alpha=.5, zorder=0)
# RT mean line on top of the band
ax_c.axvline(rt_mean, color='grey', lw=1.2, zorder=1)

# higher-EV
ax_c.plot(grand.time_ms, grand['mean'],
          color='royalblue', lw=2)
ax_c.fill_between(grand.time_ms,
                  grand['mean']-grand['sem'],
                  grand['mean']+grand['sem'],
                  color='royalblue', alpha=.25, label='_nolegend_')

# lower-EV
ax_c.plot(grand.time_ms, grand['lo_mean'],
          color='royalblue', ls='--', lw=2)
ax_c.fill_between(grand.time_ms,
                  grand['lo_mean']-grand['lo_sem'],
                  grand['lo_mean']+grand['lo_sem'],
                  color='royalblue', alpha=.20, label='_nolegend_')

ax_c.axhline(.5, ls='--', color='k')
ax_c.axvline(0,  ls='--', color='k')
ax_c.set(xlim=(0, rt_mean+100), ylim=(0, 1),
         xlabel='Time (ms)', ylabel='Fixation ratio',
         title='Choice (stimulus-locked)')

ax_c.legend(handles=line_handles, title='Gaze on', **legend_kw)
fig_choice.tight_layout()
plot_path= os.path.join(fig_dir, "gaze_choice_EE_EV.png")
fig_choice.savefig(plot_path, dpi=300, bbox_inches='tight')
fig_choice.show()

#-----------------------------------------------------------------------------------------------------------------

# plots per OV und AbsVD level

# ---------------------------------------------------------------
OV_level = {1: "low",
            2: "medium",
            3: "high"}  

block_cols = ['wheat', 'olive', 'black']           #'turquoise', 'blue', 'navy'

for OV_val, col in zip(range(1, 4), block_cols):
    label = OV_level[OV_val]
    blk_raw = raw[raw['VD_2'] == OV_val] 
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
             title=f'AbsVD-level {label}: Choice')

    # tidy legends (same style for every plot)
    # ── tidy legends (same style for every plot) ─────────────────────────────
    solid_hi  = Line2D([], [], color=col, lw=2,      label='higher EV')
    dashed_lo = Line2D([], [], color=col, ls='--', lw=2, label='lower EV')

    ax_c.legend(handles=[solid_hi, dashed_lo], frameon=False, fontsize=7)

    fig_c.savefig(os.path.join(fig_dir, f'AbsVD{label}_EE_EV.png'),dpi=300, bbox_inches='tight')




# 1) -------------------------------------------------------------
all_curves = defaultdict(list)          # keys → 'choice' | 'feed' | 'outcf'


for OV_val, col in zip(range(1, 4), block_cols):      # 1, 2, 3  (= low-medium-high)
    label   = OV_level[OV_val]
    blk_raw = raw[raw['VD_2'] == OV_val]              # ← FIX 1
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
    
# ---------------------------------------------------------------
# 2)  helper to draw one overlay figure
# ---------------------------------------------------------------
def overlay_plot(key, title, ylab, fname):
    """
    key        : 'choice', 'EE'
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
        second   = 'lo_mean' if key == 'choice' else None
        second_s = 'lo_sem'  if key == 'choice' else None
        ax.plot(df.time_ms, df[second],
                color=col, ls='--', lw=2, label=f'{label} – {dash_lbl}')
        ax.fill_between(df.time_ms,
                        df[second]-df[second_s],
                        df[second]+df[second_s],
                        color=col, alpha=.20)
        
    ax.axhline(.5, ls='--', color='k')
    ax.axvline(0 , ls='--', color='k')
    ax.set(xlim=(0, 1000), ylim=(0,1),
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
             'Choice – gaze on higher vs lower EV (AbsVD levels)',
             'Fixation ratio',
             'ALL_AbsVD_levels_EE_EV.png')

