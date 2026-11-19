from __future__ import annotations
import time, numpy as np
from numba import njit
from .parameters import Parameters
from .forces import forces_numba
from .rng import normal_pair

@njit(cache=True, fastmath=False)
def advance_numba(x, start_step, nsteps, seed, dt, gamma, noise, k_s,r0,kappa,fp,eps,sigma,rc,policy,recenter_every,offset):
    for q in range(nsteps):
        step=start_step+q
        f,_=forces_numba(x,k_s,r0,kappa,fp,eps,sigma,rc,policy)
        for i in range(x.shape[0]):
            z0,z1=normal_pair(seed,step,i)
            x[i,0]+=dt*f[i,0]/gamma+noise*z0; x[i,1]+=dt*f[i,1]/gamma+noise*z1
        if recenter_every>0 and (step+1)%recenter_every==0:
            cx=0.0; cy=0.0
            for i in range(x.shape[0]): cx+=x[i,0]; cy+=x[i,1]
            cx/=x.shape[0]; cy/=x.shape[0]; offset[0]+=cx; offset[1]+=cy
            for i in range(x.shape[0]): x[i,0]-=cx; x[i,1]-=cy
    return x

def straight(p, jitter=0.0):
    x=np.zeros((p.n_beads,2),np.float64); x[:,0]=(np.arange(p.n_beads)-p.n_bonds/2)*p.r0
    if jitter: x[:,1]=jitter*np.sin(np.linspace(0,np.pi,p.n_beads))
    return x

def _offset(recenter_every, offset):
    # Recentering is opt-in and never silent: the subtracted COM is accumulated so x+offset is the lab frame.
    if recenter_every>0 and offset is None: raise ValueError('recenter_every>0 requires an offset accumulator (np.zeros(2))')
    return np.zeros(2,np.float64) if offset is None else offset

def advance(x, step, nsteps, p, recenter_every=0, offset=None):
    return advance_numba(x,int(step),int(nsteps),int(p.seed),p.dt,p.gamma,p.noise_scale,p.k_s,p.r0,p.kappa,p.f_p,p.epsilon,p.sigma,p.cutoff,p.policy_code,int(recenter_every),_offset(recenter_every,offset))

def benchmark(p, steps=20000):
    results={}
    for name,fn in (("bruteforce",advance),("verlet",advance_verlet)):
        x=straight(p,1e-6);fn(x,0,1,p);t=time.perf_counter();fn(x,1,steps,p);elapsed=time.perf_counter()-t
        results[name]={"seconds":elapsed,"steps_per_second":steps/elapsed,"projected_hours_per_tau":p.steps_for_tau()/max(steps/elapsed,1)/3600}
    best=max(results,key=lambda k:results[k]["steps_per_second"])
    return {"backend":"numba-cpu","steps":steps,"selected":best,"variants":results,"steps_per_second":results[best]["steps_per_second"]}

@njit(cache=True, fastmath=False)
def _build_pairs(x, rc_skin, policy, pi, pj):
    n=x.shape[0];c=0;rr=rc_skin*rc_skin
    for i in range(n-1):
        for j in range(i+1,n):
            sep=j-i
            if (policy==1 and sep==1) or (policy==2 and sep<=2) or (policy==3 and sep<=3):continue
            dx=x[i,0]-x[j,0];dy=x[i,1]-x[j,1]
            if dx*dx+dy*dy<rr:pi[c]=i;pj[c]=j;c+=1
    return c

@njit(cache=True, fastmath=False)
def _forces_pairs(x,k_s,r0,kappa,fp,eps,sigma,rc,pi,pj,npairs):
    n=x.shape[0];f=np.zeros_like(x);tiny=1e-24
    for i in range(n-1):
        dx=x[i+1,0]-x[i,0];dy=x[i+1,1]-x[i,1];r=np.sqrt(dx*dx+dy*dy+tiny)
        mag=k_s*(r-r0)/r;fx=mag*dx;fy=mag*dy
        f[i,0]+=fx;f[i,1]+=fy;f[i+1,0]-=fx;f[i+1,1]-=fy
        ax=.5*fp*dx;ay=.5*fp*dy;f[i,0]+=ax;f[i,1]+=ay;f[i+1,0]+=ax;f[i+1,1]+=ay
    for i in range(n-2):
        qx=2*x[i+1,0]-x[i,0]-x[i+2,0];qy=2*x[i+1,1]-x[i,1]-x[i+2,1]
        hx=.5*kappa*qx;hy=.5*kappa*qy;f[i,0]+=hx;f[i,1]+=hy;f[i+1,0]-=2*hx;f[i+1,1]-=2*hy;f[i+2,0]+=hx;f[i+2,1]+=hy
    rc2=rc*rc;sig2=sigma*sigma
    for k in range(npairs):
        i=pi[k];j=pj[k];dx=x[i,0]-x[j,0];dy=x[i,1]-x[j,1];r2=dx*dx+dy*dy
        if r2<rc2 and r2>tiny:
            sr2=sig2/r2;sr6=sr2*sr2*sr2;sr12=sr6*sr6;c=24*eps*(2*sr12-sr6)/r2;fx=c*dx;fy=c*dy
            f[i,0]+=fx;f[i,1]+=fy;f[j,0]-=fx;f[j,1]-=fy
    return f

@njit(cache=True, fastmath=False)
def advance_verlet_numba(x,start_step,nsteps,seed,dt,gamma,noise,k_s,r0,kappa,fp,eps,sigma,rc,policy,recenter_every,skin,offset):
    n=x.shape[0];maxpairs=n*(n-1)//2;pi=np.empty(maxpairs,np.int32);pj=np.empty(maxpairs,np.int32);ref=x.copy();npairs=_build_pairs(x,rc+skin,policy,pi,pj);threshold2=(.5*skin)**2
    for q in range(nsteps):
        rebuild=False
        for i in range(n):
            dx=x[i,0]-ref[i,0];dy=x[i,1]-ref[i,1]
            if dx*dx+dy*dy>threshold2:rebuild=True;break
        if rebuild:
            npairs=_build_pairs(x,rc+skin,policy,pi,pj);ref[:,:]=x
        step=start_step+q;f=_forces_pairs(x,k_s,r0,kappa,fp,eps,sigma,rc,pi,pj,npairs)
        for i in range(n):
            z0,z1=normal_pair(seed,step,i);x[i,0]+=dt*f[i,0]/gamma+noise*z0;x[i,1]+=dt*f[i,1]/gamma+noise*z1
        if recenter_every>0 and (step+1)%recenter_every==0:
            cx=0.;cy=0.
            for i in range(n):cx+=x[i,0];cy+=x[i,1]
            cx/=n;cy/=n;offset[0]+=cx;offset[1]+=cy
            for i in range(n):x[i,0]-=cx;x[i,1]-=cy;ref[i,0]-=cx;ref[i,1]-=cy
    return x

def advance_verlet(x,step,nsteps,p,recenter_every=0,skin=0.4,offset=None):
    return advance_verlet_numba(x,int(step),int(nsteps),int(p.seed),p.dt,p.gamma,p.noise_scale,p.k_s,p.r0,p.kappa,p.f_p,p.epsilon,p.sigma,p.cutoff,p.policy_code,int(recenter_every),float(skin),_offset(recenter_every,offset))
