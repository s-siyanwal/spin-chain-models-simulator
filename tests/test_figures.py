import inspect,h5py,numpy as np,pytest
from spinchain.parameters import Parameters
from spinchain.observables import (unwrap_next,relative_bead_msd,unit_slope_fit,structure_decay_time,abs_spiral_histogram,
    standardized_kurtosis,dynamic_structure_factor,multi_origin_msd)
from spinchain.forces import total_forces
from spinchain.engine import straight
from spinchain.simulator import run
from spinchain.analysis import analyze,paper_bead_to_code
from spinchain.paper import make
from spinchain.cli import parser

def test_unwrap_next_counts_three_turns():
    th=np.linspace(0,6*np.pi,601);a=None
    for t in th:a=unwrap_next(a,(np.cos(t),np.sin(t)))
    assert abs(a-6*np.pi)<1e-9

def test_railway_msd_minimum_at_contour_transit():
    # Ideal railway motion along +x, code convention (higher index leads): bead ref-i reaches bead ref's start at lag i*r0/v.
    v=0.25;t=np.arange(400.);x=np.zeros((400,30,2));x[:,:,0]=np.arange(30)[None,:]+v*t[:,None]
    lags=np.arange(0,100);m=relative_bead_msd(x,20,20-3,lags)
    assert lags[np.argmin(m)]==12

def test_unit_slope_fit_recovers_rotational_diffusion():
    d=0.05;th=np.cumsum(np.sqrt(2*d)*np.random.default_rng(4).normal(size=200000));lags,msr=multi_origin_msd(th[:,None],2000)
    assert abs(unit_slope_fit(lags*1.0,msr)['D']/d-1)<0.1

def test_unit_slope_fit_refuses_ballistic_rotation():
    t=np.geomspace(1,1e3,60)
    assert unit_slope_fit(t,(0.3*t)**2) is None

def test_structure_decay_time_of_exponential():
    t=np.linspace(0,40,40001)
    assert abs(structure_decay_time(t,np.exp(-t/2))-2)<0.01

def test_dynamic_structure_keeps_translation_phase():
    # Rigid configuration translating at velocity u: each direction contributes |rho_0(qhat)|^2 cos(q.u t).
    # Recentering every frame would remove the phase and give S(q,t)/S(q,0)=1.
    u=np.array([0.3,0.]);base=np.c_[np.arange(5.),np.zeros(5)];x=base[None]+np.arange(50)[:,None,None]*u
    q=1.0;lags=[0,3,7];z=dynamic_structure_factor(x,[q],lags)[0];ang=np.linspace(0,2*np.pi,16,endpoint=False)
    qv=q*np.c_[np.cos(ang),np.sin(ang)];w=np.abs(np.exp(1j*base@qv.T).sum(axis=0))**2
    ref=np.array([np.sum(w*np.cos(qv@u*l))/np.sum(w) for l in lags])
    assert np.allclose(z/z[0],ref,atol=1e-12)

def test_symmetrized_kurtosis_of_constant_abs_s_is_one():
    assert abs(standardized_kurtosis(np.r_[np.full(50,2.5),np.full(50,-2.5)])-1)<1e-12

def test_abs_spiral_histogram_bins():
    assert abs_spiral_histogram([0.1,-0.1,2.6,2.4],0.25)['counts']==[2,0,0,0,0,0,0,0,0,1,1,0]

def test_paper_bead_one_is_code_leading_tip():
    q=Parameters(n_bonds=12,pe=500,xi_over_l=.2,dt=1e-5);x=straight(q);f,_=total_forces(x,q)
    assert int(np.argmax(x@f.sum(axis=0)))==paper_bead_to_code(1,q.n_beads)

def test_policy_defaults_agree():
    cli={a.dest:a.default for a in parser()._subparsers._group_actions[0].choices['headless']._actions}['policy']
    assert Parameters().nonbonded_policy==inspect.signature(make).parameters['policy'].default==cli=='lammps_default'

def _traj(tmp_path,steps=4000,every=100):
    q=Parameters(n_bonds=10,pe=200,xi_over_l=.2,dt=5e-5,seed=3);path=tmp_path/'t.h5';run(q,steps,path,sample_every=every);return q,path

def test_analysis_lag_reported_in_time_units(tmp_path):
    q,path=_traj(tmp_path)
    assert analyze(path)['com_msd']['lag_time'][0]==pytest.approx(100*q.dt)

def test_analysis_window_discards_transient(tmp_path):
    q,path=_traj(tmp_path);start=2000*q.dt/q.tau
    assert analyze(path,window_tau=(start,None),window_reason='test')['window']['samples']==21

def test_analysis_refuses_recentered_schema(tmp_path):
    path=tmp_path/'old.h5'
    with h5py.File(path,'w') as f:f.attrs['complete']=True;f.attrs['schema']=2
    with pytest.raises(ValueError):analyze(path)

def test_longest_run_above_scores_spike_as_zero_duration():
    from spinchain.observables import longest_run_above
    t=np.arange(100.);s=np.zeros(100);s[10]=5;s[40:71]=-2.5
    assert longest_run_above(t,s,2)=={"duration":30.0,"first":40,"last":70,"samples":31}

def test_paper_bead_mapping_end_points_and_railway_bead():
    assert [paper_bead_to_code(k,101) for k in (1,25,101)]==[100,76,0]  # REPRO_NOTE.md: paper bead 25 = code bead 76

def test_film_frame_selection_stride():
    from spinchain.visualize import select_frames
    assert select_frames(10,200).tolist()==list(range(10)) and select_frames(401,200).tolist()==list(range(0,401,3))
