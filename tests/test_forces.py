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
