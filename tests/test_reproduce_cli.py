import h5py,numpy as np
from spinchain.cli import main,reproduce_sample_every
from spinchain.paper import make

T1,T2=2.5e-7,5e-7  # 1263 and 2525 steps of the polymer preset: integrity sizes, not 5 tau

def _repro(d,tau,*extra):
    return main(['paper-reproduce','--regime','polymer','--tau',str(tau),'--replicates','1','--sample-every','200','--output-dir',str(d),'--i-know',*extra])

def test_default_sample_interval_is_quarter_contour():
    p=make('strong_spiral');assert reproduce_sample_every(p)==round(0.25*p.contour_time/p.dt)==1010000

def test_sample_interval_accepts_tau_units():
    p=make('strong_spiral');assert reproduce_sample_every(p,'0.01tau')==50500000

def test_sample_interval_accepts_plain_steps():
    assert reproduce_sample_every(make('strong_spiral'),'4321')==4321

def test_resume_matches_uninterrupted(tmp_path):
    a,b=tmp_path/'a',tmp_path/'b';_repro(a,T2);_repro(b,T1);_repro(b,T2,'--resume');done=make('polymer').steps_for_tau(T1)
    with h5py.File(a/'polymer-seed12001.h5','r') as f,h5py.File(b/f'polymer-seed12001.from{done}.h5','r') as g:
        assert np.array_equal(f['positions'][-1],g['positions'][-1])

def test_resume_keeps_earlier_segment(tmp_path):
    _repro(tmp_path,T1);_repro(tmp_path,T2,'--resume')
    with h5py.File(tmp_path/'polymer-seed12001.h5','r') as f:
        assert int(f['step'][-1])==make('polymer').steps_for_tau(T1)

# --- stitching paper-reproduce segments ---
import shutil,pytest
from spinchain.stitch import read_segments,write_stitched
T3=7.5e-7

def _segments(d):
    _repro(d,T1);_repro(d,T2,'--resume');_repro(d,T3,'--resume');q=make('polymer')
    return [d/'polymer-seed12001.h5',d/f'polymer-seed12001.from{q.steps_for_tau(T1)}.h5',d/f'polymer-seed12001.from{q.steps_for_tau(T2)}.h5']

def test_stitched_segments_equal_uninterrupted(tmp_path):
    segs=_segments(tmp_path/'s');_repro(tmp_path/'u',T3);out=tmp_path/'st.h5';write_stitched(segs[::-1],out)
    # The segments end at the off-grid --tau targets 1263 and 2525, so the stitched file also holds those two boundary
    # frames; every frame of the uninterrupted run must appear in it, bitwise equal.
    with h5py.File(out,'r') as f,h5py.File(tmp_path/'u'/'polymer-seed12001.h5','r') as g:
        m=np.isin(f['step'][:],g['step'][:])
        assert all(np.array_equal(f[k][:][m],g[k][:]) for k in ('step','positions','com','ee_angle'))

def test_stitch_refuses_gap(tmp_path):
    segs=_segments(tmp_path)
    with pytest.raises(ValueError,match='gap'):read_segments([segs[0],segs[2]])

def test_stitch_refuses_overlap(tmp_path):
    segs=_segments(tmp_path);dup=tmp_path/'dup.h5';shutil.copy(segs[0],dup)
    with pytest.raises(ValueError,match='overlap'):read_segments([segs[0],dup])

def test_stitch_refuses_second_origin(tmp_path):
    segs=_segments(tmp_path)
    with h5py.File(segs[1],'r+') as f:f['positions'][...]=f['positions'][...]+np.array([1.0,0.0]);f['com'][...]=f['com'][...]+np.array([1.0,0.0])
    with pytest.raises(ValueError,match='origin'):read_segments(segs[:2])

def test_seed_flag_sets_trajectory_seed(tmp_path):
    _repro(tmp_path,T1,'--seed','12013')
    with h5py.File(tmp_path/'polymer-seed12013.h5','r') as f:
        assert __import__('json').loads(f.attrs['parameters_json'])['seed']==12013
