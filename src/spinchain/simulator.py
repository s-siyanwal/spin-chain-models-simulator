from __future__ import annotations
import json,os,time,numpy as np
from .engine import straight,advance_verlet
from .io import H5Trajectory,strict_json
from .observables import end_to_end,unwrap_next
from .checkpoint import save as save_checkpoint,load as load_checkpoint

# Steps between end-to-end-angle unwrap updates (0.05 time units at dt=5e-5). The grid is global in step,
# so chunk size and resume points do not change the stored angle. Unwrapping is valid while |increment|<pi;
# the largest increment seen is written to the trajectory as ee_angle_max_increment_rad.
ANGLE_EVERY=1000

def _angle_path(checkpoint):return str(checkpoint)+'.angle.json'
def _load_angle(checkpoint,p,step):
    try:
        with open(_angle_path(checkpoint),encoding='utf8') as f:z=json.load(f)
    except (OSError,ValueError):return None
    return float(z['ee_angle']) if z.get('step')==step and z.get('fingerprint')==p.fingerprint() else None

def run(p,steps,out,sample_every=5000,chunk=5000,checkpoint=None,resume=False,overwrite=False,wall_limit=0.0,save_positions=True,angle_every=ANGLE_EVERY):
    # Checked before any file is opened: a non-positive cadence used to divide by zero, loop forever (chunk=0), or
    # "complete" with frames saved off the recorded cadence (negative sample_every).
    if int(steps)<0:raise ValueError(f'steps must be >= 0, got {steps}')
    for name,v in (('sample_every',sample_every),('chunk',chunk),('angle_every',angle_every)):
        if int(v)<1:raise ValueError(f'{name} must be a positive number of steps, got {v}')
    angle=None
    if resume:
        if not checkpoint or not os.path.exists(checkpoint):raise FileNotFoundError('resume requested but checkpoint is missing')
        x,step,wall_used=load_checkpoint(checkpoint,p);angle=_load_angle(checkpoint,p,step)
    else:x,step,wall_used=straight(p,1e-8),0,0.0
    branch='continuous' if angle is not None or not resume else 'restarted_at_resume'
    if angle is None:angle=unwrap_next(None,end_to_end(x))
    target=step+int(steps); log=H5Trajectory(out,p,overwrite=overwrite,save_positions=save_positions,sample_every=sample_every,angle_every=angle_every,start_step=step,ee_angle_branch=branch)
    t0=time.perf_counter();ok=False;max_inc=0.0
    try:
        log.append(step,p,x,angle); next_sample=((step//sample_every)+1)*sample_every
        next_progress=(int(step*p.dt/p.contour_time*2)+1)/2 if p.pe>0 else float('inf')
        while step<target:
            if wall_limit and wall_used+time.perf_counter()-t0>=wall_limit:break
            next_angle=(step//angle_every+1)*angle_every
            n=min(chunk,target-step,next_sample-step if next_sample>step else chunk,next_angle-step);advance_verlet(x,step,n,p,recenter_every=0);step+=n
            sample=step>=next_sample or step==target
            if sample or step==next_angle:
                new=unwrap_next(angle,end_to_end(x));max_inc=max(max_inc,abs(new-angle));angle=new
            if sample:
                log.append(step,p,x,angle);next_sample+=sample_every
                if checkpoint:
                    log.f.flush()  # frames up to this checkpoint are on disk before the checkpoint names them
                    save_checkpoint(checkpoint,p,x,step,wall_used+time.perf_counter()-t0)
                    strict_json(_angle_path(checkpoint),{"step":step,"fingerprint":p.fingerprint(),"ee_angle":angle})
            contours=step*p.dt/p.contour_time if p.pe>0 else 0
            if contours+1e-12>=next_progress:
                print(f'progress: step={step:,} contours={contours:.3f} tau={step*p.dt/p.tau:.6g}',flush=True);next_progress+=0.5
        ok=step>=target
    finally:log.close(complete=ok,ee_angle_max_increment_rad=max_inc)
    return {"completed":ok,"step":step,"target":target,"wall_seconds":wall_used+time.perf_counter()-t0,"output":out}
