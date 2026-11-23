import numpy as np
from spinchain.parameters import Parameters
from spinchain.forces import total_forces,conservative_energy

def test_gradient():
    p=Parameters(n_bonds=8,pe=0,xi_over_l=.5,nonbonded_policy='paper_all_pairs',dt=1e-5)
    rng=np.random.default_rng(3);x=np.c_[np.arange(9)*1.03,.2*rng.normal(size=9)]
    f,_=total_forces(x,p);num=np.zeros_like(x);h=1e-6
    for i in range(x.shape[0]):
        for a in range(2):
            xp=x.copy();xm=x.copy();xp[i,a]+=h;xm[i,a]-=h
            num[i,a]=-(conservative_energy(xp,p)-conservative_energy(xm,p))/(2*h)
    assert np.linalg.norm(f-num)/np.linalg.norm(num)<2e-5

def test_newton_and_active_total():
    p=Parameters(n_bonds=6,pe=300,xi_over_l=.5,nonbonded_policy='exclude_12_13',dt=1e-5)
    x=np.c_[np.arange(7,dtype=float),np.zeros(7)];f,_=total_forces(x,p)
    assert np.allclose(f.sum(axis=0),p.f_p*(x[-1]-x[0]))

def test_topology_policies_differ():
    x=np.array([[0.,0.],[1.,0.],[2.,0.],[3.,0.]])
    a=Parameters(n_bonds=3,pe=0,xi_over_l=.2,nonbonded_policy='paper_all_pairs',dt=1e-5)
    b=Parameters(n_bonds=3,pe=0,xi_over_l=.2,nonbonded_policy='lammps_default',dt=1e-5)
    assert total_forces(x,a)[1][2]>0 and total_forces(x,b)[1][2]==0

# --- exact properties of the force kernel ---
import pytest
from spinchain.forces import forces_numba

def test_forces_vanish_for_a_straight_chain_at_rest_length():
    p=Parameters(n_bonds=10,pe=0,xi_over_l=.5,dt=1e-5)  # lammps_default: no WCA pairs closer than 4 along the chain
    f,e=total_forces(np.c_[np.arange(11.),np.zeros(11)],p);assert np.allclose(f,0) and np.allclose(e,0)

def test_conservative_forces_have_no_net_force_or_torque():
    p=Parameters(n_bonds=12,pe=0,xi_over_l=.3,nonbonded_policy='paper_all_pairs',dt=1e-5)
    x=np.c_[np.arange(13)*0.98,0.3*np.random.default_rng(4).normal(size=13)];f,_=total_forces(x,p)
    assert np.allclose(f.sum(axis=0),0,atol=1e-9) and abs(np.sum(x[:,0]*f[:,1]-x[:,1]*f[:,0]))<1e-8

def test_wca_pair_values_at_sigma_and_at_the_cutoff():
    rc=2**(1/6)
    def pair(r):return forces_numba(np.array([[0.,0.],[r,0.]]),0.0,1.0,0.0,0.0,1.0,1.0,rc,0)  # WCA only
    f,e=pair(1.0);assert e[2]==pytest.approx(1.0) and f[0,0]==pytest.approx(-24.0) and f[1,0]==pytest.approx(24.0)
    for r in (rc*(1-1e-7),rc,1.5):  # energy and force go continuously to zero at 2^(1/6) sigma
        f,e=pair(r);assert e[2]==pytest.approx(0.0,abs=1e-12) and np.allclose(f,0,atol=1e-5)

def test_each_exclusion_policy_counts_exactly_its_pairs():
    # Every listed pair sits at r = sigma, where the WCA energy is exactly epsilon, so the energy counts pairs.
    square=np.array([[0.,0.],[1.,0.],[1.,1.],[0.,1.]])                     # separations 1 (x3) and 3 (x1) at r = 1
    rhombus=np.array([[0.,0.],[1.,0.],[.5,np.sqrt(3)/2],[1.5,np.sqrt(3)/2]])  # separations 1 (x3) and 2 (x2) at r = 1
    for pol,(a,b) in {'paper_all_pairs':(4,5),'exclude_12':(1,2),'exclude_12_13':(1,0),'lammps_default':(0,0)}.items():
        p=Parameters(n_bonds=3,pe=0,xi_over_l=.2,nonbonded_policy=pol,dt=1e-5)
        assert total_forces(square,p)[1][2]==pytest.approx(a) and total_forces(rhombus,p)[1][2]==pytest.approx(b)

def test_active_force_does_not_enter_the_conservative_energy():
    a=Parameters(n_bonds=6,pe=300,xi_over_l=.5,dt=1e-5);b=Parameters(n_bonds=6,pe=0,xi_over_l=.5,dt=1e-5)
    x=np.c_[np.arange(7)*1.02,0.1*np.sin(np.arange(7.))]
    assert conservative_energy(x,a)==conservative_energy(x,b) and np.array_equal(total_forces(x,a)[1],total_forces(x,b)[1])
