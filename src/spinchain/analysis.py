from __future__ import annotations
import json,h5py,numpy as np
from .observables import (standardized_kurtosis,multi_origin_msd,relative_bead_msd,unit_slope_fit,structure_decay_time,
    abs_spiral_histogram,kp_r2_over_l2,dynamic_structure_factor)

# Fixed analysis choices. They are documented here and are not tuned to any figure.
UNIT_SLOPE_TOL=0.1          # |d log MSR / d log t - 1| allowed in the D_r fit region (Fig. 8 rule: "gradient unity")
HIST_BIN=0.25               # |s| histogram bin width, in turns
RAILWAY_OFFSETS=(1,2,4,8)   # i in "bead j+i relative to bead j"; the paper does not list its i values
S_Q_DIRECTIONS=16           # directions of q averaged in eq. (28)

def paper_bead_to_code(k,n_beads):
    """Paper bead k (1-based, bead 1 leads under eq. (7) with r_ij=r_i-r_j) -> code index (code bead n_beads-1 leads)."""
    return n_beads-int(k)

def _clean(v):
    if isinstance(v,dict):return {k:_clean(x) for k,x in v.items()}
    if isinstance(v,(list,tuple)):return [_clean(x) for x in v]
    if isinstance(v,np.ndarray):return _clean(v.tolist())
    if isinstance(v,np.generic):v=v.item()
    if isinstance(v,float) and not np.isfinite(v):return None
    return v

def _open(path,window_tau,positions):
    with h5py.File(path,'r') as f:
        a=f.attrs
        if not bool(a.get('complete',False)):raise ValueError('trajectory is incomplete')
        if int(a.get('schema',0))<3:raise ValueError('schema<3: frames were recentered every 10000 steps with no stored offset; COM and S(q,t) cannot be repaired')
        pars=json.loads(a['parameters_json']);se=int(a['sample_every']);step=np.asarray(f['step'])
        on=np.flatnonzero(step%se==0)
        if len(on)>1 and np.any(np.diff(step[on])!=se):raise ValueError('samples are not uniformly spaced')
        tt=np.asarray(f['time_tau'])[on];lo,hi=window_tau;sel=on[(tt>=lo)&((tt<=hi) if hi is not None else True)]
        if len(sel)<2:raise ValueError(f'analysis window {window_tau} tau holds {len(sel)} samples')
        z={k:np.asarray(f[k])[sel] for k in ('step','time','time_tau','spiral','com','end_to_end','ee_angle','rg')}
        z['positions']=np.asarray(f['positions'][sel]) if positions and 'positions' in f else None
        z['attrs']={k:(v.item() if isinstance(v,np.generic) else v) for k,v in a.items() if k not in ('parameters_json',)}
    z['pars']=pars;z['sample_time']=se*pars['dt'];return z

def _lags(n,dense=32,points=160):
    m=max(1,n//2);return np.unique(np.r_[0,np.arange(1,min(dense,m)+1),np.geomspace(1,m,min(points,m)).astype(int)])

def analyze(path,include_structure=False,window_tau=(0.0,None),window_reason=None,passive_path=None):
    """Compute Figure 4/5/6/8/9/10 inputs from a schema-3 trajectory over an explicit window in tau.
    window_tau=(start,end): start discards the initial transient; for a run that never leaves a spiral, a human sets the
    window to the spiral state. No regime is inferred from beta2, |s|, or r_e."""
    z=_open(path,window_tau,include_structure);pars=z['pars'];L=pars['length'];dt_s=z['sample_time'];tau=pars['tau']
    n_beads=pars['n_bonds']+1;gamma_l=pars['gamma']*n_beads/L;kbt=pars['kbt'];fp=pars['f_p'];xi=pars['xi_over_l']*L
    s=z['spiral'];ee=z['end_to_end'];n=len(s);dur=float(z['time_tau'][-1]-z['time_tau'][0])
    warnings=[]
    if window_reason is None:window_reason='none given';warnings.append('no reason given for the analysis window; the initial transient is discarded only if window_tau[0]>0')
    if window_tau[0]<=0:warnings.append('window starts at t=0: the straight-start transient is included')
    out={"path":str(path),"schema":z['attrs']['schema'],"parameters":pars,"origin_convention":z['attrs']['origin_convention'],
        "leading_tip":z['attrs']['leading_tip'],"nonbonded_policy":pars['nonbonded_policy'],"nonbonded_policy_status":z['attrs']['nonbonded_policy_status'],
        "timestep":{"dt":pars['dt'],"status":z['attrs']['timestep_status']},"ee_angle_branch":z['attrs']['ee_angle_branch'],
        "ee_angle_max_increment_rad":z['attrs'].get('ee_angle_max_increment_rad'),
        "replication_status":("not a replication: analysis window %.3g tau < 5 tau"%dur) if dur<5 else
            "window >= 5 tau, still unconfirmed: needs the ESI timestep and the original special_bonds line",
        "window":{"start_tau":float(z['time_tau'][0]),"end_tau":float(z['time_tau'][-1]),"duration_tau":dur,"samples":n,"reason":window_reason},
        "sample_interval_time":dt_s,"sample_interval_tau":dt_s/tau,"warnings":warnings}
    out["spiral"]={"kurtosis_beta2_symmetrized":standardized_kurtosis(s),"median_abs_s":float(np.median(np.abs(s))),
        "abs_s_histogram":abs_spiral_histogram(s,HIST_BIN)}
    re2=float(np.mean(np.sum(ee**2,axis=1)))
    out["structure"]={"mean_re2_over_l2":re2/L**2,"kp_re2_over_l2_eq15":kp_r2_over_l2(pars['xi_over_l']),"mean_rg_over_l":float(np.mean(z['rg'])/L)}
    lags,msd=multi_origin_msd(z['com'])
    out["com_msd"]={"lag_time":lags*dt_s,"lag_tau":lags*dt_s/tau,"msd":msd,"com_first":z['com'][0],"com_last":z['com'][-1]}
    lags,msr=multi_origin_msd(z['ee_angle'][:,None])
    out["rotation"]={"lag_time":lags*dt_s,"msr":msr,"dr_fit_unit_slope":unit_slope_fit(lags*dt_s,msr,UNIT_SLOPE_TOL)}
    v_c=fp/gamma_l
    out["comparison_not_fit"]={"v_c":v_c,"D_t":kbt/(gamma_l*L),"D_r_a_eq23":v_c/xi,"D_r_p_rod_lower_bound":9/(4*tau),
        "v0_eq27":fp*np.sqrt(re2)/(gamma_l*L),"note":"eq. (24) D_r=D_r,p+D_r,a applies to the polymer regime only; these values are references, not fit targets"}
    pos=z['positions']
    if pos is not None:
        lg=_lags(n);j=max(1,pars['n_bonds']//4);ref=paper_bead_to_code(j,n_beads)
        out["railway"]={"paper_bead_j":j,"code_bead_j":ref,"lag_time":lg*dt_s,
            "msd_j_plus_i":{str(i):relative_bead_msd(pos,ref,ref-i,lg) for i in RAILWAY_OFFSETS if ref-i>=0},
            "self_msd_bead_j_this_run":relative_bead_msd(pos,ref,ref,lg)}
        if passive_path is not None:
            zp=_open(passive_path,window_tau,True);pp=zp['pars']
            if pp['n_bonds']!=pars['n_bonds'] or not np.isclose(zp['sample_time'],dt_s) or pp['pe']!=0:
                raise ValueError('passive reference must have Pe=0, the same N and the same sample interval')
            lp=lg[lg<=len(zp['positions'])//2]
            out["railway"]["passive_reference"]={"path":str(passive_path),"lag_time":lp*dt_s,"self_msd_bead_j":relative_bead_msd(zp['positions'],ref,ref,lp)}
        q=5*np.pi/L;sq=dynamic_structure_factor(pos,np.array([q]),lg,S_Q_DIRECTIONS)[0];c=sq/sq[0]
        out["dynamic_structure"]={"q":q,"qL":q*L,"directions":S_Q_DIRECTIONS,"lag_time":lg*dt_s,"s_over_s0":c,
            "tau_s_eq29":structure_decay_time(lg*dt_s,c),"tau_s_limits":"t=0 to t_max = floor(window samples/2) sample intervals; trapezoid rule on lag_time",
            "s_over_s0_at_upper_limit":float(c[-1])}
    return _clean(out)
