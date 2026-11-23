"""Statistical and exact checks of the noise and the integrator. All seeds are fixed, so every run is deterministic;
tolerances are about 4-5 standard errors of the quantity tested."""
import numpy as np,pytest
from numba import njit
from spinchain.rng import normal_pair
from spinchain.parameters import Parameters
from spinchain.engine import straight,advance,advance_verlet
from spinchain.observables import tangent_correlation
from spinchain.validation import von_mises_mean_cos,joint_stiffness

@njit(cache=False)
def _draw(seed,nsteps,nbeads):
    out=np.empty((nsteps,nbeads,2))
    for s in range(nsteps):
        for b in range(nbeads):
            z0,z1=normal_pair(seed,s,b);out[s,b,0]=z0;out[s,b,1]=z1
    return out

def test_counter_rng_is_standard_normal_and_uncorrelated():
    z=_draw(12001,20000,50);v=z.ravel();se=1/np.sqrt(v.size)  # 2e6 draws
    assert abs(v.mean())<5*se and abs(v.var()-1)<5*np.sqrt(2)*se
    assert abs(np.mean((v-v.mean())**4)/v.var()**2-3)<0.03
    for a,b in ((z[...,0],z[...,1]),(z[:-1],z[1:]),(z[:,:-1],z[:,1:])):  # pair partner, next step, next bead
        assert abs(np.corrcoef(a.ravel(),b.ravel())[0,1])<0.005
    assert normal_pair(1,5,3)==normal_pair(1,5,3) and normal_pair(1,5,3)!=normal_pair(2,5,3)

def test_centre_of_mass_diffusion_is_exact():
    # Internal forces cancel (Newton's third law), so at Pe=0 the COM step is the mean of the bead noises: each component
    # of a k-step COM increment has variance 2 kBT k dt/(gamma n_beads), whatever the chain's shape does.
    p=Parameters(n_bonds=3,pe=0.0,xi_over_l=10/3,seed=777);x=straight(p,1e-8);k=200;com=[]
    for i in range(20001):advance_verlet(x,i*k,k,p);com.append(x.mean(axis=0))
    var=np.diff(np.array(com),axis=0).var(axis=0);assert np.allclose(var/(2*p.kbt*k*p.dt/(p.gamma*p.n_beads)),1.0,atol=0.05)

def test_active_force_moves_the_centre_of_mass_by_f_p_r_e_over_gamma_n():
    # Same seed, so the same noise. Internal forces cancel in the COM, so COM(active) - COM(passive) is exactly the sum over
    # steps of dt f_p R_e/(gamma n_beads), with R_e the active chain's end-to-end vector before each step.
    a=Parameters(n_bonds=20,pe=1000,xi_over_l=2.0,seed=3);b=Parameters(**{**a.__dict__,'pe':0.0})
    xa,xb=straight(a,1e-8),straight(b,1e-8);drift=np.zeros(2)
    for step in range(1000):
        drift+=a.dt*a.f_p*(xa[-1]-xa[0])/(a.gamma*a.n_beads);advance(xa,step,1,a);advance(xb,step,1,b)
    assert np.allclose(xa.mean(0)-xb.mean(0),drift,rtol=0,atol=1e-10)
    v_c_t=1000*a.dt*a.f_p/a.gamma_l  # contour speed f_p/gamma_l times t; R_e <= L, so the drift is at most this
    assert drift[0]>0 and 0.8<drift[0]/v_c_t<1.001  # code bead n_beads-1 leads: a chain along +x moves to +x

def test_free_chain_bending_follows_the_boltzmann_distribution():
    # 3 bonds and no WCA pairs under lammps_default: two exactly independent joints, so <cos theta> = I1(K)/I0(K), K=5.
    p=Parameters(n_bonds=3,pe=0.0,xi_over_l=10/3,seed=777);x=straight(p,1e-8);c=[]
    for k in range(8200):
        advance_verlet(x,k*1000,1000,p)
        if k>=200:c.append(tangent_correlation(x)[1])
    assert joint_stiffness(p)==pytest.approx(5.0) and abs(np.mean(c)-von_mises_mean_cos(5.0))<0.008
