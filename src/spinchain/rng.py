import numpy as np
from numba import njit

@njit(cache=True)
def splitmix64(x):
    z=np.uint64(x)+np.uint64(0x9E3779B97F4A7C15)
    z=(z^(z>>np.uint64(30)))*np.uint64(0xBF58476D1CE4E5B9)
    z=(z^(z>>np.uint64(27)))*np.uint64(0x94D049BB133111EB)
    return z^(z>>np.uint64(31))

@njit(cache=True)
def uniform53(key):
    z=splitmix64(key)
    return float(z>>np.uint64(11))*(1.0/9007199254740992.0)

@njit(cache=True)
def normal_pair(seed,step,bead):
    key=np.uint64(seed)^((np.uint64(step)+np.uint64(1))*np.uint64(0xD2B74407B1CE6E93))^((np.uint64(bead)+np.uint64(1))*np.uint64(0xCA5A826395121157))
    u1=max(uniform53(key),1.0/9007199254740992.0)
    u2=uniform53(key^np.uint64(0x9E3779B97F4A7C15))
    rad=np.sqrt(-2.0*np.log(u1));ang=2.0*np.pi*u2
    return rad*np.cos(ang),rad*np.sin(ang)
