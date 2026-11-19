"""Formation gate evaluator (not collected by pytest). Reads schema-3 trajectories only; not a paper figure.

Pass test, and the only one: a run of consecutive saved frames with |s| > 2 lasting at least 2 contour lengths
(t_last - t_first of the run, contour time L*gamma_l/f_p), AND the stored COM equals the mean of the stored positions
with no frame after the start at |COM| ~ 0. A single-frame spike has duration 0. No other morphology measure is used.
Usage: python tests/formation_gate.py [--name gate] [--no-films] RUN.h5 [RUN.h5 ...]
"""
import argparse,os,sys,numpy as np
sys.path.insert(0,os.path.join(os.path.dirname(__file__),'..','src'))
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from spinchain.visualize import load,com_reset_check,visualize
from spinchain.observables import longest_run_above
from spinchain.io import strict_json

TITLE='formation gate, not a paper figure';S_MIN=2.0;MIN_CONTOURS=2.0
ap=argparse.ArgumentParser();ap.add_argument('paths',nargs='+');ap.add_argument('--name',default='gate');ap.add_argument('--no-films',action='store_true');args=ap.parse_args()
paths=args.paths;out_dir=os.path.dirname(paths[0]) or '.';rows=[]
fig,axes=plt.subplots(len(paths),1,figsize=(11,3.2*len(paths)),dpi=100,squeeze=False);fig.suptitle(TITLE)
for ax,path in zip(axes[:,0],paths):
    d=load(path);p=d['pars'];gl=p['gamma']*(p['n_bonds']+1)/p['length'];ct=p['length']*gl/p['f_p'];tc=d['time']/ct
    run=longest_run_above(tc,d['spiral'],S_MIN);check=com_reset_check(d)
    com_ok=check['com_equals_mean_of_positions'] and not check['com_reset_detected'];passed=run['duration']>=MIN_CONTOURS and com_ok
    row={"trajectory":path,"policy":p['nonbonded_policy'],"seed":p['seed'],"pe":p['pe'],"xi_over_l":p['xi_over_l'],"steps":int(d['step'][-1]),
        "saved_frames":len(d['step']),"contour_time":ct,"contours_run":float(tc[-1]),"duration_tau":float(d['time_tau'][-1]),
        "max_abs_s":float(np.max(np.abs(d['spiral']))),"longest_abs_s_gt_2":{**run,"duration_contours":run['duration'],
        "t_start_contours":None if run['first']<0 else float(tc[run['first']]),"t_end_contours":None if run['last']<0 else float(tc[run['last']])},
        "com_check":check,"gate_passed":bool(passed)}
    rows.append(row);print(row)
    ax.plot(tc,d['spiral'],color='k',lw=0.8);[ax.axhline(y,color='r',ls='--',lw=0.8) for y in (S_MIN,-S_MIN)]
    if run['first']>=0:ax.axvspan(tc[run['first']],tc[run['last']],color='orange',alpha=0.3,label=f"longest |s|>2 run: {run['duration']:.2f} contours")
    ax.set_title(f"{p['nonbonded_policy']}, seed {p['seed']}, Pe={p['pe']:g}, xi_P/L={p['xi_over_l']:g}: max|s|={row['max_abs_s']:.2f}, "
        f"gate {'PASSED' if passed else 'failed'} (needs |s|>2 for >= 2 contours)",fontsize=10)
    ax.set_xlabel('t / contour time (L gamma_l / f_p)');ax.set_ylabel('spiral number s')
    if run['first']>=0:ax.legend(loc='upper left',fontsize=8)
    if not args.no_films:visualize(path,out_dir,title=TITLE,window_reason='formation gate: whole run, transient included')
fig.tight_layout();fig.savefig(os.path.join(out_dir,args.name+'_s_t.png'));strict_json(os.path.join(out_dir,args.name+'_result.json'),rows)
print('ANY PASSED:',any(r['gate_passed'] for r in rows))
