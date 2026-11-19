from __future__ import annotations
import json, os, tempfile, h5py, numpy as np
from .parameters import Parameters

def strict_json(path,obj):
    def clean(v):
        if isinstance(v,dict):return {str(k):clean(x) for k,x in v.items()}
        if isinstance(v,(list,tuple)):return [clean(x) for x in v]
        if isinstance(v,np.generic):v=v.item()
        if isinstance(v,float) and not np.isfinite(v):return None
        return v
    path=os.fspath(path); os.makedirs(os.path.dirname(path) or '.',exist_ok=True)
    fd,tmp=tempfile.mkstemp(dir=os.path.dirname(path) or '.',suffix='.tmp'); os.close(fd)
    try:
        with open(tmp,'w',encoding='utf8') as f:json.dump(clean(obj),f,indent=2,sort_keys=True,allow_nan=False)
        os.replace(tmp,path)
    finally:
        if os.path.exists(tmp):os.unlink(tmp)

SCHEMA=3
ORIGIN_CONVENTION=("lab frame; origin = centre of the initial straight chain; positions are never recentered and there are "
    "no periodic images, so com and positions are unwrapped and every frame shares one origin")
LEADING_TIP=("code bead n_beads-1 leads: the kernel applies +f_p*(r[i+1]-r[i]) per bond. Eq. (7) with r_ij=r_i-r_j "
    "leads with paper bead 1, so paper bead k (1-based) = code index n_beads-k")
POLICY_STATUS={"lammps_default":"hypothesis: LAMMPS special_bonds lj 0 0 0 assumed; the paper's input script/ESI is unavailable",
    "paper_all_pairs":"printed eq. (5): all pairs j>i","exclude_12":"variant","exclude_12_13":"variant"}
TIMESTEP_STATUS="unpublished: dt is the code's stability guard 0.2*gamma/k_s, not the per-run ESI timestep"

class H5Trajectory:
    def __init__(self,path,p:Parameters,overwrite=False,save_positions=True,sample_every=None,angle_every=None,start_step=0,ee_angle_branch='continuous'):
        if os.path.exists(path) and not overwrite: raise FileExistsError(path)
        os.makedirs(os.path.dirname(path) or '.',exist_ok=True); self.f=h5py.File(path,'w'); self.save_positions=save_positions
        a=self.f.attrs;a['complete']=False; a['schema']=SCHEMA; a['fingerprint']=p.fingerprint(); a['parameters_json']=json.dumps(p.canonical(),sort_keys=True)
        a['origin_convention']=ORIGIN_CONVENTION;a['leading_tip']=LEADING_TIP;a['nonbonded_policy']=p.nonbonded_policy;a['nonbonded_policy_status']=POLICY_STATUS[p.nonbonded_policy]
        a['timestep_status']=TIMESTEP_STATUS;a['start_step']=int(start_step);a['ee_angle_branch']=ee_angle_branch
        if sample_every is not None:a['sample_every']=int(sample_every)
        if angle_every is not None:a['ee_angle_unwrap_every_steps']=int(angle_every)
        self.ds={}
        shapes={'step':(0,),'time':(0,),'time_tau':(0,),'spiral':(0,),'com':(0,2),'end_to_end':(0,2),'ee_angle':(0,),'rg':(0,)}
        if save_positions:shapes['positions']=(0,p.n_beads,2)
        for name,shape in shapes.items():
            tail=shape[1:]; self.ds[name]=self.f.create_dataset(name,shape=shape,maxshape=(None,)+tail,chunks=(min(256,max(1,256)),)+tail,dtype='f8' if name!='step' else 'i8',compression='gzip',compression_opts=1,shuffle=True)
    def append(self,step,p,x,ee_angle=None):
        from .observables import spiral_number,end_to_end,radius_gyration,unwrap_next
        ee=end_to_end(x);ang=unwrap_next(None,ee) if ee_angle is None else ee_angle
        row={'step':step,'time':step*p.dt,'time_tau':step*p.dt/p.tau,'spiral':spiral_number(x),'com':np.mean(x,axis=0),'end_to_end':ee,'ee_angle':ang,'rg':radius_gyration(x)}
        if self.save_positions:row['positions']=x
        for k,v in row.items(): d=self.ds[k]; d.resize(d.shape[0]+1,axis=0); d[-1]=v
    def close(self,complete=True,**attrs):
        for k,v in attrs.items():self.f.attrs[k]=v
        self.f.attrs['complete']=bool(complete); self.f.flush(); self.f.close()
