from __future__ import annotations
import math
from .parameters import Parameters
from .io import POLICY_STATUS

FIG2={"polymer":(200.0,0.2,100,12001),"weak_spiral":(1000.0,0.2,100,12002),"strong_spiral":(5000.0,0.2,100,12003)}
def make(regime,policy='lammps_default',dt=5e-5):
    pe,xi,n,seed=FIG2[regime];return Parameters(n_bonds=n,pe=pe,xi_over_l=xi,seed=seed,nonbonded_policy=policy,dt=dt)
def audit(p,tau_units=5.0,replicates=3,measured_sps=None):
    steps=p.steps_for_tau(tau_units);total=steps*replicates
    return {"regime_parameters":p.canonical(),"tau_units_per_replicate":tau_units,"replicates":replicates,"steps_per_replicate":steps,"total_steps":total,"estimated_hours":None if not measured_sps else total/measured_sps/3600,
        # No policy is reported as paper-faithful: the original special_bonds line is unknown (AUDIT.md).
        "nonbonded_policy":p.nonbonded_policy,"nonbonded_policy_status":POLICY_STATUS[p.nonbonded_policy],"limitations":["Main paper says timestep was adjusted; exact per-run timestep is in ESI and is not embedded here.","LAMMPS special_bonds defaults exclude 1-2 and 1-3 LJ and scale 1-4 LJ to zero; lammps_default is therefore the primary hypothesis, but ESI/input scripts are required for certainty.","Quantitative phase-map replication requires independent runs over the published grid, not only three Figure-2 points."]}
