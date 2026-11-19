import shutil,h5py,numpy as np,pytest
from spinchain.parameters import Parameters
from spinchain.simulator import run
from verdict import summarize,verdict

# Synthetic series: contour time 1, sample interval 0.25 contour, 400 samples (100 contours).
def _series(s):
    rng=np.random.default_rng(1);n=len(s);t=0.25*np.arange(n)
    return summarize(t,s,np.cumsum(rng.normal(size=(n,2)),axis=0),np.cumsum(rng.normal(size=n)),contour_time=1.0,tau=1e4)

def test_predominant_description():
    rng=np.random.default_rng(2);s=np.r_[0.3*rng.normal(size=80),-4+0.1*rng.normal(size=320)]
    assert _series(s)['descriptions']==['predominant']

def test_intermittent_description():
    rng=np.random.default_rng(3);s=0.3*rng.normal(size=400);s[200:212]=4+0.1*rng.normal(size=12)  # 2.75 contours
    assert _series(s)['descriptions']==['intermittent']

def test_absent_description():
    rng=np.random.default_rng(4);s=0.3*rng.normal(size=400);s[[50,150,250]]=5.0  # single-sample spikes
    assert _series(s)['descriptions']==['absent']

def test_short_run_refused(tmp_path):
    # N=10: tau = L^3 gamma_l / 4 = 275 time units, so 44000 steps at dt=5e-5 is 0.008 tau.
    q=Parameters(n_bonds=10,pe=200,xi_over_l=.2,dt=5e-5,seed=3);path=tmp_path/'short.h5';run(q,44000,path,sample_every=4000)
    with pytest.raises(ValueError,match='shorter than 5 tau'):verdict(path)

@pytest.fixture(scope='module')
def long_run(tmp_path_factory):
    # N=4: tau = 64*1.25/4 = 20 time units; 2.2e6 steps = 5.5 tau. Integrity size only.
    q=Parameters(n_bonds=4,pe=50,xi_over_l=.2,dt=5e-5,seed=9);path=tmp_path_factory.mktemp('v')/'long.h5'
    run(q,2200000,path,sample_every=20000);return path

def test_window_drops_first_half_tau(long_run):
    assert verdict(long_run)['window']['start_tau']==0.5

def test_missing_origin_convention_refused(long_run,tmp_path):
    p=tmp_path/'noorigin.h5';shutil.copy(long_run,p)
    with h5py.File(p,'r+') as f:del f.attrs['origin_convention']
    with pytest.raises(ValueError,match='origin_convention'):verdict(p)

def test_com_not_mean_of_positions_refused(long_run,tmp_path):
    p=tmp_path/'reset.h5';shutil.copy(long_run,p)
    with h5py.File(p,'r+') as f:f['com'][5]=np.zeros(2)
    with pytest.raises(ValueError,match='COM'):verdict(p)

from verdict import describe

def test_descriptions_at_the_rule_boundaries():
    assert describe(0.5,1.0,10.0)==[]                       # fraction exactly 1/2: neither predominant nor intermittent
    assert describe(0.9,2.5,10.0)==[]                       # mostly wound, but beta2 >= 2
    assert describe(0.9,None,10.0)==[]                      # undefined beta2 never gives predominant
    assert describe(0.6,1.0,1.0)==['predominant','absent']  # the rules are not mutually exclusive
    assert describe(0.2,1.0,2.0)==['absent']                # a plateau must be longer than 2 contours
