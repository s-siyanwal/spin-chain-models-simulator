import pytest
from spinchain.parameters import Parameters

def test_timestep_guard():
    with pytest.raises(ValueError):Parameters(dt=1e-3)
def test_dimensionless_mapping():
    p=Parameters(n_bonds=100,pe=5000,xi_over_l=.2);assert p.kappa==20 and p.f_p==.5 and p.flexure==25000
