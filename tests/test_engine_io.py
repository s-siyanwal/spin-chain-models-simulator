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

# --- engine agreement, integrity refusals and helpers ---
from spinchain.io import strict_json

def _hairpin(n_half=15,gap=1.05):
    """Chain folded back on itself: strands 1.05 apart, so beads 5+ apart along the chain start inside the WCA range."""
    a=np.c_[np.arange(n_half,dtype=float),np.zeros(n_half)];b=np.c_[np.arange(n_half-1,-1,-1,dtype=float),np.full(n_half,gap)]
    return np.vstack([a,b])

@pytest.mark.parametrize('policy',['lammps_default','paper_all_pairs'])
def test_verlet_list_engine_agrees_with_all_pairs_reference(policy):
    from spinchain.engine import advance_verlet
    q=Parameters(n_bonds=29,pe=1000,xi_over_l=.2,seed=11,nonbonded_policy=policy);a=_hairpin();b=a.copy()
    advance(a,0,20000,q);advance_verlet(b,0,20000,q)
    assert np.allclose(a,b,rtol=0,atol=1e-9)

def test_checkpoint_refuses_corrupt_state(tmp_path):
    q=p();path=tmp_path/'c.npz';save(path,q,np.full((q.n_beads,2),np.nan),5)
    with pytest.raises(ValueError):load(path,q)
    save(path,q,np.zeros((3,2)),5)
    with pytest.raises(ValueError):load(path,q)

def test_wall_limited_run_is_marked_incomplete_and_refused_by_analysis(tmp_path):
    path=tmp_path/'w.h5';r=run(p(),1000,path,sample_every=10,wall_limit=1e-12);assert not r['completed']
    with h5py.File(path,'r') as f:assert not bool(f.attrs['complete'])
    with pytest.raises(ValueError,match='incomplete'):analyze(path)

def test_resume_needs_a_checkpoint_of_the_same_physics(tmp_path):
    ck=tmp_path/'c.npz'
    with pytest.raises(FileNotFoundError):run(p(),10,tmp_path/'r.h5',sample_every=5,checkpoint=str(ck),resume=True)
    run(p(),20,tmp_path/'a.h5',sample_every=10,checkpoint=str(ck))
    with pytest.raises(ValueError,match='fingerprint'):
        run(Parameters(**{**p().__dict__,'pe':3.0}),20,tmp_path/'b.h5',sample_every=10,checkpoint=str(ck),resume=True)
    assert not (tmp_path/'b.h5').exists()

def test_strict_json_writes_null_for_non_finite_values(tmp_path):
    path=tmp_path/'a'/'b.json';strict_json(path,{'x':np.float64(np.nan),'y':[np.int64(3),float('inf')],'z':(1,2)})
    assert json.loads(path.read_text())=={'x':None,'y':[3,None],'z':[1,2]}

def test_analysis_records_policy_timestep_status_and_replication_status(tmp_path):
    q=p();path=tmp_path/'a.h5';run(q,40,path,sample_every=10);z=analyze(path)
    assert z['nonbonded_policy']=='lammps_default' and z['nonbonded_policy_status'].startswith('hypothesis')
    assert z['timestep']['status'].startswith('unpublished') and z['replication_status'].startswith('not a replication')
