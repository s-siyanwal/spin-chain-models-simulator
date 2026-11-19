"""Figures for the stitched 5tau trajectory (not collected by pytest). Not a paper figure; not a replication.

Writes, into --out-dir:
  result_panel.png  s(t) over the full run with the dropped first 0.5 tau shaded, the window |s| histogram (taken
                    from verdict.json so the figure and the verdict agree), and the lab-frame COM path.
  longest_plateau.gif  saved frames of the longest |s|>2 plateau in the window only, box following the wound chain.
Usage: python tests/result_figures.py STITCHED.h5 VERDICT.json [--out-dir DIR]
"""
import argparse,json,os,sys,numpy as np
sys.path.insert(0,os.path.join(os.path.dirname(__file__),'..','src'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from spinchain.visualize import load,plateau_film
from spinchain.observables import longest_run_above
from spinchain.io import strict_json

TITLE='5τ trajectory, not a paper figure'
ap=argparse.ArgumentParser();ap.add_argument('path');ap.add_argument('verdict');ap.add_argument('--out-dir',default='output/reproduction/result');a=ap.parse_args()
v=json.load(open(a.verdict,encoding='utf8'));d=load(a.path);p=d['pars'];os.makedirs(a.out_dir,exist_ok=True)
ct=v['contour_time'];tt=d['time_tau'];s=d['spiral'];com=d['com'];start=v['window']['start_tau']
win=np.flatnonzero(tt>=start);run=longest_run_above(d['time'][win]/ct,s[win],2.0)
first,last=int(win[run['first']]),int(win[run['last']])
if abs(run['duration']-v['longest_plateau_contours'])>1e-9:raise SystemExit('plateau recomputed here disagrees with verdict.json')

fig=plt.figure(figsize=(15,9),dpi=100);fig.suptitle(f"{TITLE}  |  Pe={p['pe']:g}  ξ_P/L={p['xi_over_l']:g}  N={p['n_bonds']}  seed={p['seed']}  policy {p['nonbonded_policy']} (unresolved vs lammps_default)",fontsize=11)
ax=fig.add_subplot(2,1,1);ax.plot(tt,s,color='k',lw=0.3);ax.axvspan(tt[0],start,color='0.85',label=f'dropped: first {start-tt[0]:g} τ')
ax.axvline(start,color='b',lw=1,label='window start');[ax.axhline(y,color='r',ls='--',lw=0.8) for y in (2,-2)]
ax.axvspan(tt[first],tt[last],color='orange',alpha=0.25,label=f"longest |s|>2 plateau: {run['duration']:.1f} contours")
ax.set_xlabel('t / τ');ax.set_ylabel('spiral number s');ax.set_xlim(tt[0],tt[-1]);ax.legend(loc='upper right',fontsize=8)
ax.set_title(f"s(t), {len(s):,} saved frames every {v['window']['sample_interval_contours']:g} contour;  window: |s|>2 fraction {v['fraction_abs_s_gt_2']:.3f}, "
    f"β₂ {v['beta2_symmetrized']:.3f}, description: {', '.join(v['descriptions']) or 'none'}",fontsize=10)
ax=fig.add_subplot(2,2,3);h=v['abs_s_histogram'];e=np.array(h['edges']);ax.bar(e[:-1],h['density'],width=h['bin_width'],align='edge',color='0.5',edgecolor='k',lw=0.4)
ax.set_xlabel('|s|');ax.set_ylabel('p(|s|)');ax.set_title(f"|s| histogram, window {start:g}–{v['window']['end_tau']:g} τ (bin {h['bin_width']:g})",fontsize=10)
ax=fig.add_subplot(2,2,4);ax.plot(com[:,0],com[:,1],color='red',lw=0.4);ax.plot(*com[0],'o',color='k',ms=4,label='start');ax.plot(*com[-1],'s',color='k',ms=4,label='end')
ax.set_aspect('equal',adjustable='datalim');ax.set_xlabel('x / r0 (lab frame)');ax.set_ylabel('y / r0');ax.legend(fontsize=8);ax.set_title('centre-of-mass path, full run (never recentered)',fontsize=10)
fig.tight_layout(rect=(0,0,1,0.96));fig.savefig(os.path.join(a.out_dir,'result_panel.png'));plt.close(fig)
film=plateau_film(a.path,os.path.join(a.out_dir,'longest_plateau.gif'),first,last,title=TITLE,follow=True,trail_from_first=True)
out={"panel":os.path.join(a.out_dir,'result_panel.png'),"film":film,"plateau_saved_frames":[first,last],"plateau_tau":[float(tt[first]),float(tt[last])],
    "plateau_contours":run['duration'],"plateau_samples":run['samples'],"com_displacement_full_run":float(np.linalg.norm(com[-1]-com[0]))}
strict_json(os.path.join(a.out_dir,'result_figures.json'),out);print(json.dumps(out,indent=1))
