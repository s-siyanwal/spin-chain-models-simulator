"""Integrity smoke run (not collected by pytest). NOT A REPLICATION: far shorter than 5 tau after a transient.

Runs one propelled Figure-2 preset and its passive (Pe=0) counterpart with the same N, xi_p/L, seed and dt, then
reads both trajectories back through spinchain.analysis (which only reads the HDF5 files) and prints a summary.
Usage: python tests/smoke_integrity.py --regime polymer --steps 4000 --sample-every 100
"""
import argparse,os,sys,numpy as np
sys.path.insert(0,os.path.join(os.path.dirname(__file__),'..','src'))
from spinchain.paper import make,FIG2
from spinchain.parameters import Parameters
from spinchain.simulator import run
from spinchain.analysis import analyze

a=argparse.ArgumentParser();a.add_argument('--regime',choices=FIG2,default='polymer');a.add_argument('--steps',type=int,default=4000)
a.add_argument('--sample-every',type=int,default=100);a.add_argument('--out-dir',default='output/smoke');a=a.parse_args()
cpu_s=2*a.steps/1.2e5
if cpu_s>180:print(f'refused: {2*a.steps:,} steps ~{cpu_s:.0f} s CPU exceeds the smoke budget');raise SystemExit(2)
act=make(a.regime);pas=Parameters(**{**act.__dict__,'pe':0.0})
pa=os.path.join(a.out_dir,f'{a.regime}-{a.steps}.h5');pp=os.path.join(a.out_dir,f'{a.regime}-passive-{a.steps}.h5')
for p,path in ((act,pa),(pas,pp)):print(run(p,a.steps,path,sample_every=a.sample_every,overwrite=True))
z=analyze(pa,include_structure=True,window_tau=(0.0,None),window_reason='smoke: whole run, transient not discarded',passive_path=pp)
print('\n'+z['replication_status']);print('policy:',z['nonbonded_policy'],'-',z['nonbonded_policy_status']);print('dt:',z['timestep'])
c0,c1=np.array(z['com_msd']['com_first']),np.array(z['com_msd']['com_last'])
print(f"unwrapped COM: first={c0} last={c1} |displacement|={np.linalg.norm(c1-c0):.4g} (lab frame, no recentering)")
with __import__('h5py').File(pa,'r') as f:s=np.asarray(f['spiral']);ang=np.asarray(f['ee_angle'])
print(f"spiral number s: first={s[0]:.3g} max|s|={np.max(np.abs(s)):.3g} median|s|={np.median(np.abs(s)):.3g}")
print(f"unwrapped ee angle: range={ang.min():.4g}..{ang.max():.4g} rad, max unwrap increment={z['ee_angle_max_increment_rad']:.3g} rad")
r=z['railway'];lt=np.array(r['lag_time']);ref=np.array(r['passive_reference']['self_msd_bead_j']);v=act.f_p*act.r0/act.gamma
print(f"railway: paper bead j={r['paper_bead_j']} (code {r['code_bead_j']}), bead speed f_p r0/gamma={v:.4g} r0 per time unit, max lag={lt[-1]:.4g}")
for i,m in r['msd_j_plus_i'].items():
    # lag 0 excluded: the passive self-MSD is 0 there. Reported, not classified.
    m=np.array(m)[:len(ref)];ratio=m[1:]/ref[1:];k=1+int(np.argmin(ratio))
    print(f"  i={i}: expected transit i*r0/v={int(i)*act.r0/v:.4g}; min MSD_(j+i)/MSD_passive = {ratio[k-1]:.3g} at lag {lt[k]:.4g}"
          f" (MSD_(j+i)={m[k]:.4g}, passive={ref[k]:.4g})"+("  [at window edge: transit not reached]" if k==len(ref)-1 else ""))
d=z['dynamic_structure'];print(f"S(q,t)/S(q,0) at qL={d['qL']:.4g}: value at upper limit {d['s_over_s0_at_upper_limit']:.3g}; tau_S(eq29, truncated)={d['tau_s_eq29']:.4g}")
print('D_r unit-slope fit:',z['rotation']['dr_fit_unit_slope']);print('comparison (not fit):',z['comparison_not_fit'])
