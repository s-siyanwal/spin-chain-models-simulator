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

# --- exact properties of the existing observables ---
import pytest

def _arc(n_bonds,alpha):
    phi=alpha*np.arange(n_bonds);return np.vstack([[0.,0.],np.cumsum(np.c_[np.cos(phi),np.sin(phi)],axis=0)])

def test_spiral_number_of_a_straight_chain_is_zero():
    assert spiral_number(np.c_[np.arange(10.),np.zeros(10)])==0.0

def test_spiral_number_counts_the_turns_of_an_open_polygon():
    # M unit bonds turning 2*pi/M at each of the M-1 joints wind (M-1)/M turns.
    for m in (6,12,40):assert spiral_number(_arc(m,2*np.pi/m))==pytest.approx((m-1)/m)

def test_spiral_number_symmetries():
    x=np.cumsum(np.random.default_rng(5).normal(size=(30,2)),axis=0);s=spiral_number(x)
    R=np.array([[np.cos(0.7),-np.sin(0.7)],[np.sin(0.7),np.cos(0.7)]])
    assert spiral_number(x@R.T+3.0)==pytest.approx(s,abs=1e-12)  # rotation and translation
    assert spiral_number(x*[1,-1])==pytest.approx(-s,abs=1e-12)    # the mirror image winds the other way
    assert spiral_number(x[::-1])==pytest.approx(-s,abs=1e-12)     # relabelling i -> N-i flips s (REPRO_NOTE.md)

def test_radius_of_gyration_and_bonds_of_simple_shapes():
    for n in (2,5,101):assert radius_gyration(np.c_[np.arange(n,dtype=float),np.zeros(n)])==pytest.approx(np.sqrt((n*n-1)/12))
    assert np.allclose(bonds(_arc(7,0.3)),1.0)

def test_kratky_porod_limits():
    a=1e3;assert kp_r2_over_l2(a)==pytest.approx(1-1/(3*a)+1/(12*a*a),abs=1e-9)  # rod: 1 - L/(3 xi_P) + ...
    assert kp_r2_over_l2(1e-4)==pytest.approx(2e-4,rel=1e-3)                       # flexible: 2 xi_P/L

def test_streaming_moments_without_symmetrization():
    x=np.random.default_rng(7).exponential(size=5000);m=Moments4()
    for z in x:m.update(z)
    assert m.kurtosis==pytest.approx(standardized_kurtosis(x,symmetrize=False),rel=1e-9)
    assert m.m1==pytest.approx(x.mean()) and m.m2/m.n==pytest.approx(x.var())

def test_kurtosis_undefined_for_tiny_or_constant_samples():
    m=Moments4()
    for z in (1.,2.,3.):m.update(z)
    assert np.isnan(m.kurtosis) and np.isnan(standardized_kurtosis([2.,2.,2.,2.],symmetrize=False))

def test_multi_origin_msd_of_a_random_walk():
    x=np.cumsum(np.random.default_rng(8).normal(size=(200000,2)),axis=0);lags,msd=multi_origin_msd(x,100)
    assert lags[0]==1 and lags[-1]==100 and np.all(np.diff(lags)>0)
    assert np.allclose(msd/lags,2.0,rtol=0.08)  # unit-variance steps in each of two components: MSD = 2*lag

def test_dynamic_structure_factor_normalisation():
    assert np.allclose(dynamic_structure_factor(np.zeros((5,1,2)),[0.7,2.0],[0]),1.0)    # one bead: S(q,0)=1
    assert np.allclose(dynamic_structure_factor(np.zeros((5,6,2)),[1.3],[0,2]),6.0)      # six coincident beads: S=N

def test_abs_spiral_histogram_density_integrates_to_one():
    h=abs_spiral_histogram(np.random.default_rng(9).normal(scale=2,size=1000),0.25)
    assert np.sum(np.array(h['density'])*h['bin_width'])==pytest.approx(1.0) and sum(h['counts'])==1000
    e=abs_spiral_histogram([],0.25);assert e['counts']==[0] and e['edges']==[0.0,0.25]

def test_unit_slope_fit_exact_line_and_too_few_points():
    assert unit_slope_fit([1,2,3],[2,4,6]) is None  # two segments < min_segments=3
    f=unit_slope_fit(np.arange(1.,6.),0.4*np.arange(1.,6.));assert f['D']==pytest.approx(0.2) and f['points']==5

def test_active_msd_eq25_limits():
    d_t,v0,d_r=0.01,0.5,0.002
    assert active_msd_eq25(np.array([0.1]),d_t,v0,d_r)[0]==pytest.approx(4*d_t*0.1+v0**2*0.01,rel=1e-4)   # 4 D_t t + v0^2 t^2
    assert active_msd_eq25(np.array([1e7]),d_t,v0,d_r)[0]==pytest.approx((4*d_t+2*v0**2/d_r)*1e7,rel=1e-3)  # diffusive

def test_longest_run_above_without_a_run_and_with_a_tie():
    t=np.arange(10.);assert longest_run_above(t,np.zeros(10),1.0)=={"duration":0.0,"first":-1,"last":-1,"samples":0}
    s=np.zeros(10);s[1:3]=3;s[6:8]=-3  # two runs of equal duration: the first is kept
    assert longest_run_above(t,s,1.0)=={"duration":1.0,"first":1,"last":2,"samples":2}

def test_structure_decay_time_of_a_constant_is_half_the_window():
    t=np.linspace(0,8,81);assert structure_decay_time(t,np.ones_like(t))==pytest.approx(4.0)

def test_unwrap_next_first_call_and_branch_crossings():
    assert unwrap_next(None,(1.0,0.0))==0.0
    assert unwrap_next(2*np.pi-0.1,(np.cos(0.1),np.sin(0.1)))==pytest.approx(2*np.pi+0.1)
    assert unwrap_next(-3.1,(np.cos(np.pi-0.05),np.sin(np.pi-0.05)))==pytest.approx(-np.pi-0.05)
