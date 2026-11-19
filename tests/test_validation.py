"""Analytic and synthetic-trajectory tests for the passive worm-like-chain check (validation.py) and its observables."""
import json,math,numpy as np,pytest
from spinchain.observables import turning_angles,tangent_correlation,spiral_number
from spinchain.validation import (von_mises_mean_cos,von_mises_pdf,joint_stiffness,discrete_persistence_length,
    passive_wlc_check)
from spinchain.parameters import Parameters
from spinchain.forces import forces_numba
from spinchain.paper import make

def _series_ratio(K):
    """Independent reference for I1(K)/I0(K): the power series of both functions."""
    i0=sum((K/2)**(2*m)/math.factorial(m)**2 for m in range(80))
    i1=sum((K/2)**(2*m+1)/(math.factorial(m)*math.factorial(m+1)) for m in range(80))
    return i1/i0

def _arc(n_bonds,alpha):
    """Unit-bond chain whose consecutive bonds turn by alpha (bond k points at angle k*alpha)."""
    phi=alpha*np.arange(n_bonds);return np.vstack([[0.0,0.0],np.cumsum(np.c_[np.cos(phi),np.sin(phi)],axis=0)])

def test_von_mises_mean_cos_matches_power_series_and_tables():
    for K in (1e-3,0.5,1.0,5.0,10.0,40.0):assert von_mises_mean_cos(K)==pytest.approx(_series_ratio(K),abs=1e-13)
    # standard table values: I0(1)=1.266065878, I1(1)=0.565159104, I0(10)=2815.716628, I1(10)=2670.988304
    assert von_mises_mean_cos(1.0)==pytest.approx(0.565159104/1.266065878,rel=1e-8)
    assert von_mises_mean_cos(10.0)==pytest.approx(2670.988304/2815.716628,rel=1e-8)

def test_von_mises_mean_cos_limits():
    assert von_mises_mean_cos(1e-4)==pytest.approx(0.5e-4,rel=1e-6)                     # K/2 for small K
    assert von_mises_mean_cos(500.0)==pytest.approx(1-1/1000-1/(8*500**2),abs=1e-8)     # 1-1/2K-1/8K^2 for large K

def test_von_mises_pdf_is_normalised_and_has_the_exact_mean_cos():
    t=np.linspace(-np.pi,np.pi,200001)
    for K in (0.5,10.0,80.0):
        assert np.trapezoid(von_mises_pdf(t,K),t)==pytest.approx(1.0,abs=1e-10)
        assert np.trapezoid(np.cos(t)*von_mises_pdf(t,K),t)==pytest.approx(von_mises_mean_cos(K),abs=1e-10)

def test_joint_stiffness_of_the_figure2_presets_is_ten():
    # kappa = (xi_P/L) kBT L = 0.2*100 = 20 for every Figure-2 preset, so K = kappa r0^2/(2 kBT) = 10.
    for regime in ('polymer','weak_spiral','strong_spiral'):assert joint_stiffness(make(regime))==pytest.approx(10.0)
    assert joint_stiffness(Parameters(n_bonds=20,pe=0.0,xi_over_l=1.0))==pytest.approx(10.0)

def test_kernel_joint_energy_is_k_kbt_one_minus_cos():
    # The printed bending energy of the force kernel, per joint, is K kBT (1 - cos theta) at unit bond length.
    p=Parameters(n_bonds=3,pe=0.0,xi_over_l=1.0)
    for theta in (0.0,0.3,1.1,2.5):
        _,e=forces_numba(_arc(3,theta),p.k_s,p.r0,p.kappa,0.0,p.epsilon,p.sigma,p.cutoff,p.policy_code)
        assert e[1]==pytest.approx(2*joint_stiffness(p)*p.kbt*(1-np.cos(theta)),rel=1e-12,abs=1e-14)

def test_discrete_persistence_length():
    assert discrete_persistence_length(1e4)==pytest.approx(2e4-1,abs=0.01)   # xi_P - r0 for large K
    assert discrete_persistence_length(10.0)==pytest.approx(-1/np.log(2670.988304/2815.716628),rel=1e-8)

def test_turning_angles_of_an_arc():
    x=_arc(12,0.2)
    assert np.allclose(turning_angles(x),0.2) and np.allclose(turning_angles(x[::-1]),-0.2)  # reversal flips the sign
    assert spiral_number(x)==pytest.approx(turning_angles(x).sum()/(2*np.pi))
    assert turning_angles(np.stack([x,_arc(12,-0.1)])).shape==(2,11)

def test_tangent_correlation_of_a_constant_curvature_arc_is_cos_n_alpha():
    assert np.allclose(tangent_correlation(_arc(30,0.15)),np.cos(0.15*np.arange(30)),atol=1e-12)

def test_tangent_correlation_is_invariant_under_rigid_motion_and_reversal():
    x=np.cumsum(np.random.default_rng(0).normal(size=(25,2)),axis=0);c=tangent_correlation(x)
    R=np.array([[np.cos(1.3),-np.sin(1.3)],[np.sin(1.3),np.cos(1.3)]])
    assert c[0]==pytest.approx(1.0)
    assert np.allclose(tangent_correlation(x@R.T+7.0),c,atol=1e-12) and np.allclose(tangent_correlation(x[::-1]),c,atol=1e-12)
    with pytest.raises(ValueError):tangent_correlation(x,max_sep=24)

def test_synthetic_von_mises_chains_give_u_to_the_n():
    # Synthetic trajectory: 20000 chains of 20 bonds with i.i.d. von Mises turning angles (K=10), random orientation.
    # C(n) must be u^n; the standard error at n=19 is about 0.004.
    K=10.0;rng=np.random.default_rng(1);th=rng.vonmises(0.0,K,size=(20000,19))
    phi=np.c_[np.zeros(20000),np.cumsum(th,axis=1)]+rng.uniform(0,2*np.pi,size=(20000,1))
    b=np.stack([np.cos(phi),np.sin(phi)],axis=-1);x=np.concatenate([np.zeros((20000,1,2)),np.cumsum(b,axis=1)],axis=1)
    assert np.allclose(turning_angles(x),th,atol=1e-9)
    assert np.allclose(tangent_correlation(x),von_mises_mean_cos(K)**np.arange(20),atol=0.02)

def test_passive_wlc_check_pipeline(tmp_path):
    # Integrity size only (5 bonds, 2 seeds): run -> HDF5 -> analysis reader -> statistics -> JSON and PNG.
    res,traj=tmp_path/'res',tmp_path/'traj'
    passive_wlc_check(out_dir=str(res),traj_dir=str(traj),n_bonds=5,xi_over_l=4.0,seeds=(1,2),steps=400_000,
                      burn_steps=40_000,sample_every=400)
    z=json.loads((res/'passive_wlc_check.json').read_text())
    assert (res/'passive_wlc_check.png').stat().st_size>10_000 and len(list(traj.glob('*.h5')))==2
    assert z['model']['joint_stiffness_K']==pytest.approx(10.0) and z['run']['seeds']==[1,2]
    assert z['run']['frames_per_seed_after_burn_in']==[901,901] and 'not a paper figure' in z['note']
    r=z['results'];assert r['tangent_correlation_mean'][0]==pytest.approx(1.0) and r['z_vs_exact'][0] is None
    assert abs(r['mean_cos']-z['theory']['mean_cos_exact_I1_over_I0'])<0.01

def test_passive_wlc_check_refuses_one_seed(tmp_path):
    with pytest.raises(ValueError,match='two seeds'):passive_wlc_check(out_dir=str(tmp_path),traj_dir=str(tmp_path),seeds=(1,))
