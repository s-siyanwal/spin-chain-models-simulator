"""Passive-limit (Pe = 0) check of the bending term against exact results. A code check: not a paper figure, not a
replication.

With Pe = 0 the model is a discrete worm-like chain. With bond length |b_i| = r0, the printed bending energy of one joint
(Eq. 4) is

    (kappa/4) |b_i - b_(i+1)|^2 = (kappa r0^2 / 2) (1 - cos theta_i),

so in equilibrium, for a free chain without excluded-volume contacts, the turning angles theta_i are independent with

    p(theta) = exp(K cos theta) / (2 pi I0(K)),    K = kappa r0^2 / (2 kBT).

Hence <cos theta> = I1(K)/I0(K) = u, and the exact tangent correlation of the discrete chain is C(n) = <t_i . t_(i+n)> = u^n.
The continuum worm-like chain has C(s) = exp(-s/xi_P) with xi_P = kappa/kBT (Eq. 8). The discrete chain's persistence
length -r0/ln(u) tends to xi_P - r0 for large K, so the two curves differ by a few per cent at K = 10. Bond-length
fluctuations (k_S = 4000 kBT/r0^2) change the effective K by about 1e-4, below the statistical error of this check.

The check runs the production path (simulator.run -> schema-3 HDF5 -> the analysis reader) on a chain short enough for
every bending mode to relax within the run. The 100-bond Figure-2 chain has the same joint stiffness (kappa = 20 kBT r0,
K = 10), but its slowest bending mode relaxes in about 2e4 time units (4e8 steps), far longer than a quick check allows.
"""
from __future__ import annotations
import os,platform,time
import numpy as np
from . import __version__
from .parameters import Parameters
from .simulator import run
from .analysis import _open,_clean
from .observables import turning_angles,tangent_correlation
from .io import strict_json,POLICY_STATUS,TIMESTEP_STATUS

TITLE='passive-limit check of the bending term, not a paper figure'
SEEDS=(1,2,3,4,5,6,7,8)

def _quad_grid(n=20001):return np.linspace(0.0,np.pi,n)

def von_mises_mean_cos(K):
    """Exact <cos theta> = I1(K)/I0(K) for p(theta) ~ exp(K cos theta). The trapezoid rule on [0, pi] is spectrally
    accurate for this smooth periodic integrand (agrees with the power series to ~1e-16); exp(K(cos-1)) avoids overflow."""
    t=_quad_grid();w=np.exp(float(K)*(np.cos(t)-1.0));return float(np.trapezoid(w*np.cos(t),t)/np.trapezoid(w,t))

def von_mises_pdf(theta,K):
    """p(theta) = exp(K cos theta)/(2 pi I0(K)) on (-pi, pi], evaluated without overflow."""
    t=_quad_grid();i0_scaled=np.trapezoid(np.exp(float(K)*(np.cos(t)-1.0)),t)/np.pi  # exp(-K) I0(K)
    return np.exp(float(K)*(np.cos(np.asarray(theta,float))-1.0))/(2*np.pi*i0_scaled)

def joint_stiffness(p:Parameters):
    """K = kappa r0^2/(2 kBT): the joint energy (kappa/4)|b_i-b_(i+1)|^2 is (kappa r0^2/2)(1-cos theta) at |b|=r0."""
    return p.kappa*p.r0**2/(2.0*p.kbt)

def discrete_persistence_length(K,r0=1.0):
    """-r0/ln(I1(K)/I0(K)): the decay length of C(n)=u^n. Tends to 2K r0 - r0, i.e. xi_P - r0, for large K."""
    return -float(r0)/np.log(von_mises_mean_cos(K))

def passive_wlc_check(out_dir='results/validation',traj_dir='output/validation',n_bonds=20,xi_over_l=1.0,seeds=SEEDS,
                      steps=8_000_000,burn_steps=2_000_000,sample_every=2000,policy='lammps_default',
                      name='passive_wlc_check',figure=True):
    """Run Pe=0 trajectories for each seed, drop the burn-in, and compare the turning-angle distribution and C(n) with
    the exact results above. Error bars are standard errors over the independent seeds. Writes <name>.json (and
    <name>.png if figure) into out_dir; the HDF5 trajectories go to traj_dir (git-ignored by *.h5)."""
    if len(seeds)<2:raise ValueError('need at least two seeds for a standard error over seeds')
    if not 0<=burn_steps<steps:raise ValueError('burn_steps must be in [0, steps)')
    t_start=time.perf_counter();edges=np.linspace(-np.pi,np.pi,121);per=[];paths=[]
    for seed in seeds:
        p=Parameters(n_bonds=n_bonds,pe=0.0,xi_over_l=xi_over_l,seed=int(seed),nonbonded_policy=policy)
        path=os.path.join(traj_dir,f'{name}-N{n_bonds}-seed{int(seed)}.h5');paths.append(path)
        run(p,steps,path,sample_every=sample_every,overwrite=True)
        z=_open(path,(burn_steps*p.dt/p.tau,None),True)  # the analysis reader: complete, schema 3, uniform cadence
        pos=z['positions'];dens,_=np.histogram(turning_angles(pos).ravel(),edges,density=True)
        per.append({"seed":int(seed),"fingerprint":p.fingerprint(),"frames":int(len(pos)),"first_step":int(z['step'][0]),
            "tangent_correlation":tangent_correlation(pos),"angle_density":dens})
    K=joint_stiffness(p);u=von_mises_mean_cos(K);n=np.arange(n_bonds);nseed=len(per)
    C=np.array([s['tangent_correlation'] for s in per]);H=np.array([s['angle_density'] for s in per])
    c_mean=C.mean(0);c_se=C.std(0,ddof=1)/np.sqrt(nseed);h_mean=H.mean(0);h_se=H.std(0,ddof=1)/np.sqrt(nseed)
    exact=u**n;z=np.full(n_bonds,np.nan);z[1:]=(c_mean[1:]-exact[1:])/c_se[1:]
    centres=0.5*(edges[1:]+edges[:-1])
    out={"note":"code check of the bending term in the passive limit (Pe=0); not a paper figure and not a replication",
        "title":TITLE,
        "model":{"n_bonds":n_bonds,"pe":0.0,"xi_over_l":xi_over_l,"kappa":p.kappa,"kbt":p.kbt,"r0":p.r0,
            "joint_stiffness_K":K,"dt":p.dt,"timestep_status":TIMESTEP_STATUS,"nonbonded_policy":policy,
            "nonbonded_policy_status":POLICY_STATUS[policy]},
        "run":{"seeds":[s['seed'] for s in per],"steps_per_seed":int(steps),"burn_in_steps":int(burn_steps),
            "sample_every":int(sample_every),"frames_per_seed_after_burn_in":[s['frames'] for s in per],
            "trajectories_not_in_git":paths,"fingerprints":{str(s['seed']):s['fingerprint'] for s in per}},
        "theory":{"mean_cos_exact_I1_over_I0":u,"discrete_persistence_length_r0":discrete_persistence_length(K,p.r0),
            "continuum_xi_P_r0":p.kappa/p.kbt,"tangent_correlation_exact":exact,
            "tangent_correlation_continuum":np.exp(-n*p.r0*p.kbt/p.kappa)},
        "results":{"separation_bonds":n,"tangent_correlation_mean":c_mean,"tangent_correlation_se_over_seeds":c_se,
            "z_vs_exact":z,"max_abs_z":float(np.nanmax(np.abs(z))),"mean_cos":float(c_mean[1]),"mean_cos_se":float(c_se[1]),
            "per_seed_mean_cos":[float(s['tangent_correlation'][1]) for s in per],
            "angle_bin_edges":edges,"angle_density_mean":h_mean,"angle_density_se_over_seeds":h_se,
            "angle_density_exact_at_bin_centres":von_mises_pdf(centres,K)},
        "provenance":_versions()|{"runtime_seconds":time.perf_counter()-t_start}}
    os.makedirs(out_dir,exist_ok=True)
    if figure:out["figure"]=_figure(out,os.path.join(out_dir,name+'.png'))
    out["json"]=os.path.join(out_dir,name+'.json');strict_json(out["json"],_clean(out))
    return out

def _versions():
    import numba,h5py
    return {"spinchain_version":__version__,"python":platform.python_version(),"numpy":np.__version__,
        "numba":numba.__version__,"h5py":h5py.__version__,"system":f"{platform.system()} {platform.machine()}"}

def render_figure(json_path,png_path=None):
    """Redraw the figure from a saved <name>.json without rerunning any simulation (same numbers, same plot)."""
    import json
    with open(json_path,encoding='utf8') as f:z=json.load(f)
    return _figure(z,png_path or os.path.splitext(json_path)[0]+'.png')

def _figure(z,path):
    """Two panels: turning-angle density and tangent correlation, simulation (+-1 s.e. over seeds) against exact."""
    from matplotlib.figure import Figure  # object-oriented API: no global backend change, safe inside notebooks
    from matplotlib.ticker import MaxNLocator
    m,r,th=z['model'],z['run'],z['theory'];res=z['results'];K=m['joint_stiffness_K']
    fig=Figure(figsize=(12,5.2),dpi=150,layout='constrained');ax=fig.subplots(1,2)
    e=np.asarray(res['angle_bin_edges']);c=0.5*(e[1:]+e[:-1]);h=np.asarray(res['angle_density_mean']);hs=np.asarray(res['angle_density_se_over_seeds'])
    ok=h>0;g=np.linspace(-1.7,1.7,601)
    a=ax[0];a.errorbar(c[ok],h[ok],yerr=hs[ok],fmt='o',ms=3.2,color='C0',ecolor='C0',elinewidth=0.8,capsize=0,label='simulation (mean of seeds)')
    a.plot(g,von_mises_pdf(g,K),color='k',lw=1.4,label=rf'exact: $e^{{K\cos\theta}}/2\pi I_0(K)$, $K={K:g}$')
    a.set_yscale('log');a.set_xlim(-1.7,1.7);a.set_ylim(1e-4,4);a.set_xlabel(r'turning angle $\theta$ between consecutive bonds (rad)')
    a.set_ylabel(r'probability density $p(\theta)$');a.set_title('(a) bending-angle distribution',fontsize=11);a.legend(fontsize=8.5,loc='lower center',frameon=False)
    n=np.asarray(res['separation_bonds']);cm=np.asarray(res['tangent_correlation_mean']);cs=np.asarray(res['tangent_correlation_se_over_seeds'])
    a=ax[1];a.plot(n,th['tangent_correlation_exact'],color='k',lw=1.4,label=rf"exact, discrete chain: $u^n$, $\ell_p=-r_0/\ln u={th['discrete_persistence_length_r0']:.2f}\,r_0$")
    a.plot(n,th['tangent_correlation_continuum'],color='0.45',ls='--',lw=1.2,label=rf"continuum WLC: $e^{{-n r_0/\xi_P}}$, $\xi_P=\kappa/k_BT={th['continuum_xi_P_r0']:g}\,r_0$")
    a.errorbar(n,cm,yerr=cs,fmt='o',ms=3.6,color='C3',ecolor='C3',elinewidth=0.9,capsize=2,label='simulation (mean of seeds)')
    a.set_yscale('log');a.xaxis.set_major_locator(MaxNLocator(integer=True))
    a.set_xlabel(r'bond separation $n$');a.set_ylabel(r'$C(n)=\langle \mathbf{t}_i\cdot\mathbf{t}_{i+n}\rangle$')
    a.set_title(r'(b) tangent correlation along the chain',fontsize=11);a.legend(fontsize=8.5,loc='lower left',frameon=False)
    a.text(0.97,0.97,rf"$\langle\cos\theta\rangle$ = {res['mean_cos']:.5f} $\pm$ {res['mean_cos_se']:.5f}"+"\n"+
        rf"exact $I_1(K)/I_0(K)$ = {th['mean_cos_exact_I1_over_I0']:.5f}"+"\n"+rf"max $|z|$ over $n$ = {res['max_abs_z']:.2f}",
        transform=a.transAxes,ha='right',va='top',fontsize=8.5)
    steps=r['steps_per_seed']-r['burn_in_steps']
    fig.suptitle(f"{TITLE}\nPe = 0, N = {m['n_bonds']} bonds, $\\kappa$ = {m['kappa']:g} $k_BT\\,r_0$ (K = {K:g}, the joint stiffness of the "
        f"Fig. 2 presets), policy {m['nonbonded_policy']}, dt = {m['dt']:g} (unpublished)\n{len(r['seeds'])} seeds x {steps:,} steps after a "
        f"{r['burn_in_steps']:,}-step burn-in; error bars: $\\pm$1 standard error over seeds",fontsize=9.5)
    os.makedirs(os.path.dirname(path) or '.',exist_ok=True);fig.savefig(path)
    return path
