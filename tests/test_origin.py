import h5py,numpy as np,pytest
from spinchain.parameters import Parameters
from spinchain.engine import straight,advance_verlet
from spinchain.simulator import run

# Propelled 10-bond chain. 10000 steps is the old silent-recenter boundary of advance_verlet.
def p():return Parameters(n_bonds=10,pe=200,xi_over_l=.2,dt=5e-5,nonbonded_policy='lammps_default',seed=7)

def test_stored_com_not_reset_on_recenter_boundary(tmp_path):
    path=tmp_path/'a.h5';run(p(),10000,path,sample_every=10000)
    with h5py.File(path,'r') as f:com=np.asarray(f['com'][-1])
    assert np.linalg.norm(com)>1e-3

def test_stored_positions_share_engine_origin(tmp_path):
    q=p();path=tmp_path/'a.h5';run(q,12000,path,sample_every=3000,chunk=2500)
    x=straight(q,1e-8);advance_verlet(x,0,12000,q)
    with h5py.File(path,'r') as f:last=np.asarray(f['positions'][-1])
    assert np.array_equal(last,x)

def test_stored_com_is_mean_of_stored_positions(tmp_path):
    path=tmp_path/'a.h5';run(p(),12000,path,sample_every=3000)
    with h5py.File(path,'r') as f:com=np.asarray(f['com']);pos=np.asarray(f['positions'])
    assert np.array_equal(com,pos.mean(axis=1))

def test_recenter_without_offset_refused():
    q=p();x=straight(q)
    with pytest.raises(ValueError):advance_verlet(x,0,10,q,recenter_every=5)

def test_recenter_offset_restores_lab_frame():
    q=p();a=straight(q);b=a.copy();off=np.zeros(2)
    advance_verlet(a,0,3000,q);advance_verlet(b,0,3000,q,recenter_every=1000,offset=off)
    assert np.allclose(b+off,a,rtol=0,atol=1e-9)

def test_resume_positions_match_uninterrupted(tmp_path):
    q=p();ck=tmp_path/'c.npz'
    run(q,6000,tmp_path/'full.h5',sample_every=1000)
    run(q,3000,tmp_path/'a.h5',sample_every=1000,checkpoint=str(ck));run(q,3000,tmp_path/'b.h5',sample_every=1000,checkpoint=str(ck),resume=True)
    with h5py.File(tmp_path/'full.h5','r') as f,h5py.File(tmp_path/'b.h5','r') as g:
        assert np.array_equal(f['positions'][-1],g['positions'][-1])

def test_resume_ee_angle_matches_uninterrupted(tmp_path):
    q=p();ck=tmp_path/'c.npz'
    run(q,6000,tmp_path/'full.h5',sample_every=1000)
    run(q,3000,tmp_path/'a.h5',sample_every=1000,checkpoint=str(ck));run(q,3000,tmp_path/'b.h5',sample_every=1000,checkpoint=str(ck),resume=True)
    with h5py.File(tmp_path/'full.h5','r') as f,h5py.File(tmp_path/'b.h5','r') as g:
        assert np.array_equal(f['ee_angle'][3:],g['ee_angle'][:])

def test_run_chunk_size_does_not_change_trajectory(tmp_path):
    q=p();run(q,5000,tmp_path/'a.h5',sample_every=1000,chunk=5000);run(q,5000,tmp_path/'b.h5',sample_every=1000,chunk=333)
    with h5py.File(tmp_path/'a.h5','r') as f,h5py.File(tmp_path/'b.h5','r') as g:
        assert np.array_equal(f['positions'][:],g['positions'][:])
