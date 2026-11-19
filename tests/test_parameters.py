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

import dataclasses,math

def test_derived_quantities_of_the_strong_spiral_preset():
    p=Parameters(n_bonds=100,pe=5000,xi_over_l=.2)
    assert p.gamma_l==pytest.approx(1.01) and p.tau==pytest.approx(252500) and p.contour_time==pytest.approx(202.0)
    assert p.noise_scale==pytest.approx(math.sqrt(2*5e-5)) and p.cutoff==pytest.approx(2**(1/6)) and p.k_s==4000

def test_fingerprint_is_deterministic_and_records_the_policy():
    a=Parameters(seed=1);assert a.fingerprint()==Parameters(seed=1).fingerprint()
    assert Parameters(seed=2).fingerprint()!=a.fingerprint()
    assert Parameters(seed=1,nonbonded_policy='paper_all_pairs').fingerprint()!=a.fingerprint()  # AGENTS.md rule

def test_invalid_parameters_are_refused():
    for kw in ({'n_bonds':2},{'nonbonded_policy':'all'},{'xi_over_l':0},{'kbt':0},{'gamma':-1},{'sigma_over_r0':0}):
        with pytest.raises(ValueError):Parameters(**kw)

def test_parameters_are_immutable():
    with pytest.raises(dataclasses.FrozenInstanceError):Parameters().pe=1.0
