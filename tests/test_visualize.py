import h5py,numpy as np,pytest
from PIL import Image
from spinchain.parameters import Parameters
from spinchain.simulator import run
from spinchain.visualize import load,select_frames,com_reset_check,visualize

def test_refuses_file_without_origin_convention(tmp_path):
    path=tmp_path/'old.h5'
    with h5py.File(path,'w') as f:f.attrs['schema']=2;f.create_dataset('positions',data=np.zeros((2,3,2)))
    with pytest.raises(ValueError):load(path)

def test_film_frame_cap():
    assert len(select_frames(201,200))<=200

def test_gif_frames_are_saved_frames(tmp_path):
    q=Parameters(n_bonds=10,pe=200,xi_over_l=.2,dt=5e-5,seed=5);path=tmp_path/'t.h5';run(q,2400,path,sample_every=100)
    out=visualize(path,str(tmp_path/'viz'),max_frames=10)
    assert Image.open(out['gif']).n_frames==len(select_frames(25,10))

def test_reset_check_flags_recentered_com():
    pos=np.random.default_rng(0).normal(size=(4,5,2));pos-=pos.mean(axis=1,keepdims=True)
    d={'com':pos.mean(axis=1),'positions':pos,'step':np.arange(4)*10000}
    assert com_reset_check(d)['com_reset_detected']

def test_plateau_film_uses_only_requested_saved_frames(tmp_path):
    from spinchain.visualize import plateau_film
    q=Parameters(n_bonds=10,pe=200,xi_over_l=.2,dt=5e-5,seed=5);path=tmp_path/'t.h5';run(q,2400,path,sample_every=100)
    out=plateau_film(path,str(tmp_path/'p.gif'),5,17)
    assert Image.open(out['gif']).n_frames==13

def test_plateau_film_follow_mode_uses_same_saved_frames(tmp_path):
    from spinchain.visualize import plateau_film
    q=Parameters(n_bonds=10,pe=200,xi_over_l=.2,dt=5e-5,seed=5);path=tmp_path/'t.h5';run(q,2400,path,sample_every=100)
    out=plateau_film(path,str(tmp_path/'f.gif'),5,17,follow=True)
    assert Image.open(out['gif']).n_frames==13
