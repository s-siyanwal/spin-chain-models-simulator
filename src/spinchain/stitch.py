"""Stitch paper-reproduce segments (<name>.h5, <name>.from<step>.h5, ...) into one schema-3 trajectory.

Consecutive segments must meet at exactly one shared frame: a resumed segment starts at the checkpoint step, which is
also the last frame of the previous segment. That frame must agree bitwise (positions and COM), which is the one-origin
check. A later start is a gap and an earlier start is an overlap; both are refused, not repaired. Segments must share
the physics fingerprint, origin convention, sample cadence, and policy. A non-final segment may be incomplete (the
process was stopped after its last checkpoint); the stitched file is complete only if the last segment is.
"""
from __future__ import annotations
import json,os,h5py,numpy as np
from .io import SCHEMA

FIELDS=('step','time','time_tau','spiral','com','end_to_end','ee_angle','rg','positions')
SAME=('schema','fingerprint','parameters_json','origin_convention','leading_tip','sample_every','nonbonded_policy')

def _read(path):
    with h5py.File(path,'r') as f:
        a={k:(v.item() if isinstance(v,np.generic) else v) for k,v in f.attrs.items()}
        if int(a.get('schema',0))!=SCHEMA or 'origin_convention' not in a:raise ValueError(f'{path}: not a schema-{SCHEMA} lab-frame trajectory')
        d={k:np.asarray(f[k]) for k in FIELDS if k in f}
    if len(d['step'])==0:raise ValueError(f'{path}: empty segment')
    return a,d,str(path)

def read_segments(paths):
    """Return (attrs, data) for the concatenation of the given segments in step order."""
    segs=sorted((_read(p) for p in paths),key=lambda s:int(s[1]['step'][0]))
    a0,d0,_=segs[0];keys=[k for k in FIELDS if all(k in s[1] for s in segs)]
    for key in SAME:
        if any(s[0].get(key)!=a0.get(key) for s in segs):raise ValueError(f'segments disagree on {key}')
    out={k:[d0[k]] for k in keys};prev_a,prev_d,prev_p=segs[0];max_inc=float(a0.get('ee_angle_max_increment_rad',0.0))
    for a,d,p in segs[1:]:
        last=int(prev_d['step'][-1]);first=int(d['step'][0])
        if first<last:raise ValueError(f'overlap: {p} starts at step {first} but {prev_p} ends at step {last}')
        if first>last:raise ValueError(f'gap: {prev_p} ends at step {last} but {p} starts at step {first}')
        if int(a.get('start_step',first))!=first:raise ValueError(f'{p}: start_step attribute disagrees with its first frame')
        for k in ('positions','com'):
            if k in keys and not np.array_equal(prev_d[k][-1],d[k][0]):raise ValueError(f'origin: {k} at shared step {first} differs between {prev_p} and {p}')
        ang=d['ee_angle'];shift=float(prev_d['ee_angle'][-1]-ang[0])
        if a.get('ee_angle_branch')=='continuous':
            if shift!=0:raise ValueError(f'{p}: continuous ee_angle branch disagrees at shared step {first}')
        else:
            turns=shift/(2*np.pi)
            if abs(turns-round(turns))>1e-9:raise ValueError(f'{p}: restarted ee_angle branch is not a 2*pi shift')
            d={**d,'ee_angle':ang+2*np.pi*round(turns)}
        for k in keys:out[k].append(d[k][1:])
        max_inc=max(max_inc,float(a.get('ee_angle_max_increment_rad',0.0)));prev_a,prev_d,prev_p=a,d,p
    data={k:np.concatenate(v) for k,v in out.items()}
    if np.any(np.diff(data['step'])<=0):raise ValueError('stitched steps are not strictly increasing')
    attrs={**a0,'complete':bool(segs[-1][0].get('complete',False)),'start_step':int(data['step'][0]),'ee_angle_branch':'continuous',
        'ee_angle_max_increment_rad':max_inc,'stitched_from':json.dumps([s[2] for s in segs])}
    return attrs,data

def write_stitched(paths,out,overwrite=False):
    """Write the stitched trajectory to one schema-3 file that analyze() and visualize() read unchanged."""
    if os.path.exists(out) and not overwrite:raise FileExistsError(out)
    attrs,data=read_segments(paths);os.makedirs(os.path.dirname(out) or '.',exist_ok=True)
    with h5py.File(out,'w') as f:
        for k,v in data.items():f.create_dataset(k,data=v,chunks=(min(256,len(v)),)+v.shape[1:],compression='gzip',compression_opts=1,shuffle=True)
        for k,v in attrs.items():f.attrs[k]=v
    return {"output":str(out),"segments":json.loads(attrs['stitched_from']),"frames":int(len(data['step'])),
        "first_step":int(data['step'][0]),"last_step":int(data['step'][-1]),"complete":attrs['complete']}
