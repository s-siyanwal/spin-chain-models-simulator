from __future__ import annotations
import numpy as np

def spiral_number(x):
    b=np.diff(np.asarray(x),axis=0); phi=np.arctan2(b[:,1],b[:,0]); d=(np.diff(phi)+np.pi)%(2*np.pi)-np.pi
    return float(d.sum()/(2*np.pi))
def end_to_end(x): return np.asarray(x[-1]-x[0],float)
def radius_gyration(x):
    y=np.asarray(x)-np.mean(x,axis=0); return float(np.sqrt(np.mean(np.sum(y*y,axis=1))))
def bonds(x): return np.linalg.norm(np.diff(x,axis=0),axis=1)
def kp_r2_over_l2(xi_over_l):
    a=float(xi_over_l); return 2*a-2*a*a*(1-np.exp(-1/a))
def standardized_kurtosis(values,symmetrize=True):
    v=np.asarray(values,float)
    if symmetrize: v=np.concatenate((np.abs(v),-np.abs(v)))
    if v.size<4 or np.var(v)==0:return float("nan")
    d=v-v.mean(); return float(np.mean(d**4)/np.mean(d**2)**2)
def multi_origin_msd(v,max_lag=None):
    a=np.asarray(v,float); n=len(a); m=min(max_lag or n//2,n-1); lags=np.unique(np.geomspace(1,m,min(160,m)).astype(int)); out=[]
    for lag in lags: out.append(np.mean(np.sum((a[lag:]-a[:-lag])**2,axis=-1)))
    return lags,np.asarray(out)
def unwrapped_angle(v): return np.unwrap(np.arctan2(np.asarray(v)[:,1],np.asarray(v)[:,0]))
def multi_origin_msr(v,max_lag=None):
    a=unwrapped_angle(v)[:,None]; return multi_origin_msd(a,max_lag)
def dynamic_structure_factor(positions,q_values,lags,directions=16):
    x=np.asarray(positions,float); qs=np.asarray(q_values,float); result=np.empty((len(qs),len(lags)))
    angles=np.linspace(0,2*np.pi,directions,endpoint=False)
    for iq,q in enumerate(qs):
        acc=np.zeros(len(lags)); count=np.zeros(len(lags))
        for ang in angles:
            qv=q*np.array([np.cos(ang),np.sin(ang)]); rho=np.exp(1j*np.tensordot(x,qv,axes=([2],[0]))).sum(axis=1)
            for k,lag in enumerate(lags):
                if lag==0: z=rho*np.conj(rho)
                else:z=rho[:-lag]*np.conj(rho[lag:])
                acc[k]+=z.real.mean();count[k]+=1
        result[iq]=acc/count/(x.shape[1])
    return result

# --- Figure 6/8/9/10 inputs. Lags passed in are sample indices; callers convert to time units. ---
def unwrap_next(prev,v):
    """Continue an unwrapped angle: atan2(v) on the branch nearest prev. Valid while |increment| < pi."""
    a=float(np.arctan2(v[1],v[0]))
    return a if prev is None else float(prev+(a-prev+np.pi)%(2*np.pi)-np.pi)
def relative_bead_msd(positions,ref,other,lags):
    """<|r_other(t0+lag)-r_ref(t0)|^2> over origins t0 (Fig. 6). other==ref gives the self MSD of bead ref."""
    x=np.asarray(positions,float);out=[]
    for lag in lags:
        d=x[lag:,other]-x[:len(x)-lag,ref];out.append(np.mean(np.sum(d*d,axis=-1)))
    return np.asarray(out)
def unit_slope_fit(t,y,tol=0.1,min_segments=3):
    """Fit y=2*D*t on the longest contiguous run of log-log segments with |slope-1|<=tol (eq. 16/Fig. 8 rule).
    Returns None when no run of min_segments consecutive unit-slope segments exists."""
    t=np.asarray(t,float);y=np.asarray(y,float);ok=(t>0)&(y>0)&np.isfinite(y);t=t[ok];y=y[ok]
    if len(t)<min_segments+1:return None
    slope=np.diff(np.log(y))/np.diff(np.log(t));good=np.abs(slope-1)<=tol;best=(0,0);start=None
    for k,g in enumerate(np.r_[good,False]):
        if g and start is None:start=k
        if not g and start is not None:
            if k-start>best[1]-best[0]:best=(start,k)
            start=None
    if best[1]-best[0]<min_segments:return None
    a,b=best[0],best[1]+1;d=float(np.exp(np.mean(np.log(y[a:b])-np.log(2*t[a:b]))))
    return {"D":d,"t_min":float(t[a]),"t_max":float(t[b-1]),"points":int(b-a),"slope_tolerance":tol}
def structure_decay_time(t,c):
    """Eq. (29): tau_S = int t S dt / int S dt, trapezoid rule over the supplied grid; limits are t[0]..t[-1]."""
    t=np.asarray(t,float);c=np.asarray(c,float);return float(np.trapezoid(t*c,t)/np.trapezoid(c,t))
def abs_spiral_histogram(s,bin_width=0.25):
    """Histogram of |s| with fixed bins [0, bin_width, ...]; density integrates to 1."""
    a=np.abs(np.asarray(s,float));top=bin_width*max(1,int(np.ceil(a.max()/bin_width))+1) if a.size else bin_width
    edges=np.arange(0,top+bin_width/2,bin_width);counts,_=np.histogram(a,edges)
    return {"bin_width":bin_width,"edges":edges.tolist(),"counts":counts.tolist(),"density":(counts/max(a.size,1)/bin_width).tolist()}
def active_msd_eq25(t,d_t,v0,d_r):
    """Eq. (25): 4 D_t t + (2 v0^2/D_r^2)(D_r t + exp(-D_r t) - 1). Comparison curve only."""
    t=np.asarray(t,float);return 4*d_t*t+2*v0*v0/d_r**2*(d_r*t+np.exp(-d_r*t)-1)
def longest_run_above(t,values,threshold):
    """Longest run of consecutive saved samples with |value| > threshold. Duration is t[last]-t[first] of that run,
    so a single-sample spike has duration 0. Behaviour between samples is not observed."""
    t=np.asarray(t,float);a=np.abs(np.asarray(values,float))>threshold;best=(0.0,-1,-1);start=None
    for k,g in enumerate(np.r_[a,False]):
        if g and start is None:start=k
        if not g and start is not None:
            dur=float(t[k-1]-t[start])
            if dur>best[0] or best[1]<0:best=(dur,start,k-1)
            start=None
    return {"duration":best[0],"first":int(best[1]),"last":int(best[2]),"samples":0 if best[1]<0 else int(best[2]-best[1]+1)}

# --- Passive worm-like-chain checks (validation.py). positions: (n_beads,2) or (frames,n_beads,2). ---
def turning_angles(positions):
    """Signed turning angle between consecutive bonds, in (-pi, pi]: the same wrapped difference of bond angles that
    spiral_number sums. Shape (..., n_bonds-1)."""
    b=np.diff(np.asarray(positions,float),axis=-2);phi=np.arctan2(b[...,1],b[...,0])
    return (np.diff(phi,axis=-1)+np.pi)%(2*np.pi)-np.pi
def tangent_correlation(positions,max_sep=None):
    """C(n)=<t_i . t_(i+n)> for n=0..max_sep: unit bond tangents, averaged over bond pairs i and over frames. C(0)=1."""
    b=np.diff(np.asarray(positions,float),axis=-2);t=b/np.linalg.norm(b,axis=-1,keepdims=True)
    if t.ndim==2:t=t[None]
    nb=t.shape[1];m=nb-1 if max_sep is None else int(max_sep)
    if not 0<=m<nb:raise ValueError(f'max_sep must be in 0..{nb-1}, got {max_sep}')
    return np.array([float(np.mean(np.sum(t[:,:nb-n]*t[:,n:],axis=-1))) for n in range(m+1)])
