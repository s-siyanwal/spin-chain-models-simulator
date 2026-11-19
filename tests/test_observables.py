import numpy as np
from spinchain.observables import *
from spinchain.moments import Moments4

def test_spiral_number():
    t=np.linspace(0,6*np.pi,500);r=.3*t;x=np.c_[r*np.cos(t),r*np.sin(t)]
    assert abs(abs(spiral_number(x))-3.23)<.05

def test_kurtosis_gaussian():
    x=np.random.default_rng(1).normal(size=200000);assert abs(standardized_kurtosis(x,False)-3)<.05

def test_streaming_moments():
    x=np.random.default_rng(2).normal(size=10000);m=Moments4()
    for z in x:m.add_symmetrized(z)
    assert abs(m.kurtosis-standardized_kurtosis(x,True))<1e-10

def test_kp(): assert abs(kp_r2_over_l2(.2)-0.3205390358)<1e-9

def test_msd_ballistic():
    t=np.arange(500.);x=np.c_[2*t,np.zeros_like(t)];l,y=multi_origin_msd(x)
    assert abs(np.polyfit(np.log(l),np.log(y),1)[0]-2)<1e-12

def test_dynamic_structure_static():
    x=np.zeros((20,4,2));x[:]=np.c_[np.arange(4),np.zeros(4)];z=dynamic_structure_factor(x,[1.0],[0,1,5])
    assert np.allclose(z,z[:,0,None])
