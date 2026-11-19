import pytest
from spinchain.parameters import Parameters

def test_timestep_guard():
    with pytest.raises(ValueError):Parameters(dt=1e-3)
def test_dimensionless_mapping():
    p=Parameters(n_bonds=100,pe=5000,xi_over_l=.2);assert p.kappa==20 and p.f_p==.5 and p.flexure==25000

def test_r0_other_than_one_is_refused():
    # The printed bending energy kappa/4 sum|b_i-b_(i+1)|^2 is written in units r0=1: with r0=2 a simulated chain is
    # r0^3 = 8 times stiffer than xi_P/L says (measured l_p/xi_P = 7.99), so other r0 are refused rather than mislabelled.
    with pytest.raises(ValueError,match='r0'):Parameters(r0=2.0)
    assert Parameters(r0=1.0).length==100
