import json,h5py,numpy as np,pytest
from spinchain.parameters import Parameters
from spinchain.engine import straight,advance
from spinchain.simulator import run
from spinchain.analysis import analyze
from spinchain.checkpoint import save,load

def p():return Parameters(n_bonds=10,pe=2,xi_over_l=.2,dt=1e-5,nonbonded_policy='lammps_default')
def test_chunk_invariance():
    q=p();a=straight(q);b=a.copy();advance(a,0,25,q);advance(b,0,10,q);advance(b,10,15,q);assert np.array_equal(a,b)
def test_checkpoint(tmp_path):
    q=p();x=straight(q);path=tmp_path/'c.npz';save(path,q,x,17,2.5);y,s,w=load(path,q);assert np.array_equal(x,y) and s==17 and w==2.5
    q2=Parameters(**{**q.__dict__,'pe':3});
    with pytest.raises(ValueError):load(path,q2)
def test_h5_exact_and_analyze(tmp_path):
    q=p();path=tmp_path/'a.h5';r=run(q,23,path,sample_every=10,chunk=7);assert r['completed']
    with h5py.File(path,'r') as f:assert bool(f.attrs['complete']) and list(f['step'])==[0,10,20,23]
    z=analyze(path);json.dumps(z,allow_nan=False)
def test_overwrite_refusal(tmp_path):
    q=p();path=tmp_path/'a.h5';run(q,2,path,sample_every=1)
    with pytest.raises(FileExistsError):run(q,2,path,sample_every=1)
def test_verlet_chunk_invariance():
    from spinchain.engine import advance_verlet
    q=p();a=straight(q);b=a.copy();advance_verlet(a,0,25,q);advance_verlet(b,0,10,q);advance_verlet(b,10,15,q);assert np.array_equal(a,b)

@pytest.mark.parametrize('kw',[{'sample_every':0},{'sample_every':-5},{'chunk':0},{'angle_every':0},{'steps':-1}])
def test_run_refuses_non_positive_cadence_before_creating_output(tmp_path,kw):
    # sample_every=-5 used to "complete" while saving every chunk under a recorded cadence of -5, sample_every=0 raised
    # ZeroDivisionError after the output file was created, and chunk=0 looped forever without advancing.
    path=tmp_path/'bad.h5';args={'steps':10,'sample_every':5,**kw};steps=args.pop('steps')
    with pytest.raises(ValueError):run(p(),steps,path,**args)
    assert not path.exists()
