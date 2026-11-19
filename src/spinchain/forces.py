from __future__ import annotations
import numpy as np
from numba import njit
from .parameters import Parameters

@njit(cache=True, fastmath=False)
def forces_numba(x, k_s, r0, kappa, fp, eps, sigma, rc, policy):
    n=x.shape[0]; f=np.zeros_like(x); e=np.zeros(3, np.float64)
    tiny=1e-24
    for i in range(n-1):
        dx=x[i+1,0]-x[i,0]; dy=x[i+1,1]-x[i,1]; r=(dx*dx+dy*dy+tiny)**0.5
        mag=k_s*(r-r0)/r; fx=mag*dx; fy=mag*dy
        f[i,0]+=fx; f[i,1]+=fy; f[i+1,0]-=fx; f[i+1,1]-=fy
        e[0]+=0.5*k_s*(r-r0)*(r-r0)
        ax=0.5*fp*dx; ay=0.5*fp*dy
        f[i,0]+=ax; f[i,1]+=ay; f[i+1,0]+=ax; f[i+1,1]+=ay
    for i in range(n-2):
        qx=2*x[i+1,0]-x[i,0]-x[i+2,0]; qy=2*x[i+1,1]-x[i,1]-x[i+2,1]
        hx=0.5*kappa*qx; hy=0.5*kappa*qy
        f[i,0]+=hx; f[i,1]+=hy; f[i+1,0]-=2*hx; f[i+1,1]-=2*hy; f[i+2,0]+=hx; f[i+2,1]+=hy
        e[1]+=0.25*kappa*(qx*qx+qy*qy)
    rc2=rc*rc; sig2=sigma*sigma
    for i in range(n-1):
        for j in range(i+1,n):
            sep=j-i
            if (policy==1 and sep==1) or (policy==2 and sep<=2) or (policy==3 and sep<=3): continue
            dx=x[i,0]-x[j,0]; dy=x[i,1]-x[j,1]; r2=dx*dx+dy*dy
            if r2<rc2 and r2>tiny:
                sr2=sig2/r2; sr6=sr2*sr2*sr2; sr12=sr6*sr6
                c=24*eps*(2*sr12-sr6)/r2; fx=c*dx; fy=c*dy
                f[i,0]+=fx; f[i,1]+=fy; f[j,0]-=fx; f[j,1]-=fy
                e[2]+=4*eps*(sr12-sr6)+eps
    return f,e

def total_forces(x: np.ndarray, p: Parameters):
    return forces_numba(np.asarray(x,np.float64),p.k_s,p.r0,p.kappa,p.f_p,p.epsilon,p.sigma,p.cutoff,p.policy_code)

def conservative_energy(x: np.ndarray,p: Parameters):
    _,e=forces_numba(np.asarray(x,np.float64),p.k_s,p.r0,p.kappa,0.0,p.epsilon,p.sigma,p.cutoff,p.policy_code)
    return float(e.sum())
