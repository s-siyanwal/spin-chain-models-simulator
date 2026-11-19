"""Verdict on one stitched schema-3 trajectory (not collected by pytest). Not a replication.

Reads one file, refuses it unless it has origin_convention, a stored COM equal to the mean of the stored positions,
and a duration of at least 5 tau; drops the first 0.5 tau (the only window choice, recorded in the output); and writes
verdict.json with these quantities and nothing else:
  duration in tau; fraction of window samples with |s|>2; longest |s|>2 plateau in contours; symmetrized beta2;
  |s| histogram; COM MSD log-log slope; whether a unit-slope MSR region exists (D_r null if not).
The descriptions 'predominant' / 'intermittent' / 'absent' are text applied by fixed rules. They are not a pass bit
and not a claim that any figure was reproduced. No morphology cut selects samples.
Usage: python tests/verdict.py STITCHED.h5 [--out verdict.json]
"""
import argparse,json,os,sys,h5py,numpy as np
sys.path.insert(0,os.path.join(os.path.dirname(__file__),'..','src'))
from spinchain.observables import (standardized_kurtosis,abs_spiral_histogram,multi_origin_msd,unit_slope_fit,longest_run_above)
from spinchain.analysis import UNIT_SLOPE_TOL,HIST_BIN
from spinchain.io import strict_json

MIN_DURATION_TAU=5.0;DROP_TAU=0.5;S_THRESHOLD=2.0;PLATEAU_CONTOURS=2.0;HALF=0.5;BETA2_MAX=2.0
RULES={"predominant":"fraction of window samples with |s|>2 is above 1/2 and symmetrized beta2 < 2",
    "intermittent":"a |s|>2 plateau longer than 2 contours exists and the |s|>2 fraction is below 1/2",
    "absent":"no |s|>2 plateau longer than 2 contours exists"}

def describe(frac,beta2,plateau_contours):
    """All descriptions whose rule holds, in RULES order. Can be empty, or hold more than one."""
    out=[]
    if frac>HALF and beta2 is not None and beta2<BETA2_MAX:out.append('predominant')
    if plateau_contours>PLATEAU_CONTOURS and frac<HALF:out.append('intermittent')
    if not plateau_contours>PLATEAU_CONTOURS:out.append('absent')
    return out

def summarize(time,s,com,ee_angle,contour_time,tau):
    """Window quantities. time in simulation units, uniformly sampled; inputs are already windowed."""
    s=np.asarray(s,float);t=np.asarray(time,float);dt_s=float(t[1]-t[0])
    frac=float(np.mean(np.abs(s)>S_THRESHOLD));run=longest_run_above(t/contour_time,s,S_THRESHOLD);b2=standardized_kurtosis(s)
    lags,msd=multi_origin_msd(np.asarray(com,float));ok=msd>0
    slope=float(np.polyfit(np.log(lags[ok]*dt_s),np.log(msd[ok]),1)[0]) if ok.sum()>=2 else None
    lags,msr=multi_origin_msd(np.asarray(ee_angle,float)[:,None]);fit=unit_slope_fit(lags*dt_s,msr,UNIT_SLOPE_TOL)
    b2=None if not np.isfinite(b2) else float(b2)
    return {"fraction_abs_s_gt_2":frac,"longest_plateau_contours":float(run['duration']),"beta2_symmetrized":b2,
        "abs_s_histogram":abs_spiral_histogram(s,HIST_BIN),
        "com_msd_loglog_slope":slope,"com_msd_slope_rule":"least-squares slope of log MSD vs log lag over lags from one sample to half the window",
        "msr_unit_slope_region_exists":fit is not None,"D_r":None if fit is None else fit['D'],"msr_unit_slope_fit":fit,
        "descriptions":describe(frac,b2,float(run['duration']))}

def load_checked(path):
    with h5py.File(path,'r') as f:
        a=f.attrs
        if 'origin_convention' not in a:raise ValueError(f'{path}: refused, origin_convention attribute missing')
        if 'positions' not in f:raise ValueError(f'{path}: refused, no positions to check the COM against')
        pars=json.loads(a['parameters_json']);se=int(a['sample_every'])
        d={k:np.asarray(f[k]) for k in ('step','time','time_tau','spiral','com','ee_angle','positions')}
        info={"complete":bool(a.get('complete',False)),"stitched_from":a.get('stitched_from'),"nonbonded_policy":a.get('nonbonded_policy'),
            "nonbonded_policy_status":a.get('nonbonded_policy_status'),"timestep_status":a.get('timestep_status')}
    if not np.array_equal(d['com'],d['positions'].mean(axis=1)):raise ValueError(f'{path}: refused, stored COM is not the mean of the stored positions')
    dur=float(d['time_tau'][-1]-d['time_tau'][0])
    if dur<MIN_DURATION_TAU*(1-1e-9):raise ValueError(f'{path}: refused, duration {dur:.4g} tau is shorter than {MIN_DURATION_TAU:g} tau')
    on=d['step']%se==0
    if np.any(np.diff(d['step'][on])!=se):raise ValueError(f'{path}: refused, saved frames are not uniformly spaced')
    return pars,se,{k:v[on] for k,v in d.items()},dur,info

def verdict(path):
    pars,se,d,dur,info=load_checked(path);gl=pars['gamma']*(pars['n_bonds']+1)/pars['length'];ct=pars['length']*gl/pars['f_p']
    start=float(d['time_tau'][0]+DROP_TAU);w=d['time_tau']>=start
    out={"note":"not a replication; descriptions are not a claim that a figure was reproduced","input":str(path),**info,
        "parameters":{k:pars[k] for k in ('pe','xi_over_l','n_bonds','seed','dt','tau')},"contour_time":ct,
        "duration_tau":dur,"window":{"rule":f"drop the first {DROP_TAU} tau; no morphology cut","start_tau":start,
        "end_tau":float(d['time_tau'][-1]),"samples":int(w.sum()),"sample_interval_contours":se*pars['dt']/ct},
        "description_rules":RULES}
    out.update(summarize(d['time'][w],d['spiral'][w],d['com'][w],d['ee_angle'][w],ct,pars['tau']));return out

if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('path');ap.add_argument('--out',default='verdict.json');a=ap.parse_args()
    try:v=verdict(a.path)
    except ValueError as e:print(e,file=sys.stderr);raise SystemExit(2)
    strict_json(a.out,v);print(json.dumps({k:v[k] for k in ('duration_tau','fraction_abs_s_gt_2','longest_plateau_contours','beta2_symmetrized','com_msd_loglog_slope','msr_unit_slope_region_exists','D_r','descriptions')},indent=2))
