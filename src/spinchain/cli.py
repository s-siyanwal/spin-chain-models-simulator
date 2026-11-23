from __future__ import annotations
import argparse,json,os,sys,time
import numpy as np
from .parameters import Parameters
from .paper import make,audit,FIG2
from .engine import benchmark
from .simulator import run
from .analysis import analyze
from .io import strict_json
from .hardware import probe
from .taichi_backend import probe_vulkan

REPRODUCE_SAMPLE_DEFAULT='0.25contour'
def reproduce_sample_every(p,spec=None):
    """paper-reproduce cadence in steps. spec: '<x>contour' (contour time L*gamma_l/f_p), '<x>tau', or integer steps.
    Default 0.25 contour. Checkpoints are written at every saved frame."""
    spec=str(spec or REPRODUCE_SAMPLE_DEFAULT).strip().lower()
    for unit,scale in (('contours',p.contour_time),('contour',p.contour_time),('tau',p.tau)):
        if spec.endswith(unit):return max(1,int(round(float(spec[:-len(unit)])*scale/p.dt)))
    return int(spec)
def pp(x):print(json.dumps(x,indent=2,sort_keys=True,allow_nan=False,default=str))
def parser():
    p=argparse.ArgumentParser(prog='spinchain');s=p.add_subparsers(dest='cmd',required=True)
    q=s.add_parser('benchmark');q.add_argument('--steps',type=int,default=20000);q.add_argument('--regime',choices=FIG2,default='polymer');q.add_argument('--policy',choices=['lammps_default','paper_all_pairs','exclude_12','exclude_12_13'],default='lammps_default')
    q=s.add_parser('paper-audit');q.add_argument('--regime',choices=FIG2,default='polymer');q.add_argument('--tau',type=float,default=5);q.add_argument('--replicates',type=int,default=3);q.add_argument('--policy',choices=['lammps_default','paper_all_pairs','exclude_12','exclude_12_13'],default='lammps_default');q.add_argument('--benchmark-steps',type=int,default=20000)
    q=s.add_parser('headless');q.add_argument('--regime',choices=FIG2,default='polymer');q.add_argument('--steps',type=int,required=True);q.add_argument('--sample-every',type=int,default=5000);q.add_argument('--output',required=True);q.add_argument('--checkpoint');q.add_argument('--resume',action='store_true');q.add_argument('--overwrite',action='store_true');q.add_argument('--policy',choices=['lammps_default','paper_all_pairs','exclude_12','exclude_12_13'],default='lammps_default');q.add_argument('--wall-limit',type=float,default=0);q.add_argument('--no-positions',action='store_true')
    q=s.add_parser('paper-reproduce');q.add_argument('--regime',choices=FIG2,required=True);q.add_argument('--tau',type=float,default=5);q.add_argument('--replicates',type=int,default=3);q.add_argument('--output-dir',default='output/reproduction');q.add_argument('--policy',choices=['lammps_default','paper_all_pairs','exclude_12','exclude_12_13'],default='lammps_default');q.add_argument('--i-know',action='store_true');q.add_argument('--overwrite',action='store_true')
    q.add_argument('--sample-every',help="saved-frame and checkpoint interval: '<x>contour', '<x>tau', or integer steps (default 0.25contour)");q.add_argument('--resume',action='store_true',help='continue from the checkpoint into a new segment file');q.add_argument('--seed',type=int,help='base seed; replicate r uses seed+r (default: the Fig. 2 preset seed)')
    q=s.add_parser('analyze');q.add_argument('path');q.add_argument('--structure',action='store_true');q.add_argument('--json-out')
    q.add_argument('--window-start-tau',type=float,default=0.0);q.add_argument('--window-end-tau',type=float);q.add_argument('--window-reason')
    q=s.add_parser('visualize');q.add_argument('path');q.add_argument('--passive');q.add_argument('--out-dir',default='output/viz');q.add_argument('--fps',type=int,default=12);q.add_argument('--max-frames',type=int,default=200)
    q=s.add_parser('wlc-check',help='passive-limit (Pe=0) check of the bending term against exact results; not a paper figure')
    q.add_argument('--out-dir',default='results/validation');q.add_argument('--traj-dir',default='output/validation')
    s.add_parser('hardware');s.add_parser('vulkan-probe');return p

def main(argv=None):
    a=parser().parse_args(argv)
    if a.cmd=='benchmark':pp(benchmark(make(a.regime,a.policy),a.steps));return 0
    if a.cmd=='paper-audit':
        p=make(a.regime,a.policy);b=benchmark(p,a.benchmark_steps);pp({"benchmark":b,"audit":audit(p,a.tau,a.replicates,b['steps_per_second'])});return 0
    if a.cmd=='headless':pp(run(make(a.regime,a.policy),a.steps,a.output,a.sample_every,checkpoint=a.checkpoint,resume=a.resume,overwrite=a.overwrite,wall_limit=a.wall_limit,save_positions=not a.no_positions));return 0
    if a.cmd=='paper-reproduce':
        if not a.i_know:print('Refused: full reproduction requires --i-know after paper-audit.',file=sys.stderr);return 2
        code=0
        for r in range(a.replicates):
            p0=make(a.regime,a.policy);p=Parameters(**{**p0.__dict__,"seed":(p0.seed if a.seed is None else a.seed)+r});steps=p.steps_for_tau(a.tau);path=os.path.join(a.output_dir,f'{a.regime}-seed{p.seed}.h5');ck=path+'.checkpoint.npz'
            se=reproduce_sample_every(p,a.sample_every);done=0
            if a.resume:
                # A resume never reopens the earlier segment: it writes <name>.from<step>.h5 and runs only the remaining steps.
                with np.load(ck) as z:done=int(z['step'])
                path=os.path.join(a.output_dir,f'{a.regime}-seed{p.seed}.from{done}.h5')
            pp({"plan":{"total_steps":steps,"done_steps":done,"remaining_steps":steps-done,"sample_every":se,"sample_interval_tau":se*p.dt/p.tau,"sample_interval_contours":se*p.dt/p.contour_time,"output":path,"checkpoint":ck,"policy":p.nonbonded_policy,"seed":p.seed}})
            res=run(p,steps-done,path,sample_every=se,checkpoint=ck,resume=a.resume,overwrite=a.overwrite);pp(res);code=max(code,0 if res['completed'] else 2)
        return code
    if a.cmd=='analyze':
        z=analyze(a.path,a.structure,window_tau=(a.window_start_tau,a.window_end_tau),window_reason=a.window_reason);pp(z)
        if a.json_out:strict_json(a.json_out,z)
        return 0
    if a.cmd=='visualize':
        from .visualize import visualize
        pp(visualize(a.path,a.out_dir,a.passive,a.fps,a.max_frames));return 0
    if a.cmd=='wlc-check':
        # Fixed protocol (validation.passive_wlc_check defaults): reruns exactly the committed check, about a minute of CPU.
        from .validation import passive_wlc_check
        z=passive_wlc_check(out_dir=a.out_dir,traj_dir=a.traj_dir);r=z['results']
        pp({"json":z['json'],"figure":z.get('figure'),"mean_cos":r['mean_cos'],"mean_cos_se":r['mean_cos_se'],
            "mean_cos_exact":z['theory']['mean_cos_exact_I1_over_I0'],"max_abs_z":r['max_abs_z'],
            "runtime_seconds":z['provenance']['runtime_seconds'],"note":z['note']});return 0
    if a.cmd=='hardware':pp(probe());return 0
    if a.cmd=='vulkan-probe':pp(probe_vulkan());return 0
