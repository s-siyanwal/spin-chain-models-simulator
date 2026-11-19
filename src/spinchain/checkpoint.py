from __future__ import annotations
import os,tempfile,numpy as np
from .parameters import Parameters
SCHEMA=1
def save(path,p:Parameters,x,step,wall_seconds=0.0):
    os.makedirs(os.path.dirname(path) or '.',exist_ok=True); fd,tmp=tempfile.mkstemp(dir=os.path.dirname(path) or '.',suffix='.npz');os.close(fd)
    try:
        np.savez_compressed(tmp,x=np.asarray(x,np.float64),step=np.int64(step),wall_seconds=np.float64(wall_seconds),schema=np.int64(SCHEMA),fingerprint=np.array(p.fingerprint()))
        os.replace(tmp,path)
    finally:
        if os.path.exists(tmp):os.unlink(tmp)
def load(path,p:Parameters):
    with np.load(path,allow_pickle=False) as z:
        if int(z['schema'])!=SCHEMA:raise ValueError('checkpoint schema mismatch')
        if str(z['fingerprint'])!=p.fingerprint():raise ValueError('checkpoint physics fingerprint mismatch')
        x=np.array(z['x'],np.float64);step=int(z['step']);wall=float(z['wall_seconds'])
    if x.shape!=(p.n_beads,2) or not np.isfinite(x).all() or step<0:raise ValueError('invalid checkpoint state')
    return x,step,wall
