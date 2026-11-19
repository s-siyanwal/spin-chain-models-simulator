from __future__ import annotations
from dataclasses import dataclass, asdict
import hashlib, json, math

POLICIES = {"paper_all_pairs": 0, "exclude_12": 1, "exclude_12_13": 2, "lammps_default": 3}

@dataclass(frozen=True)
class Parameters:
    n_bonds: int = 100
    r0: float = 1.0
    sigma_over_r0: float = 1.0
    pe: float = 200.0
    xi_over_l: float = 0.2
    kbt: float = 1.0
    gamma: float = 1.0
    epsilon: float = 1.0
    spring_factor: float = 4000.0
    dt: float = 5.0e-5
    nonbonded_policy: str = "lammps_default"  # hypothesis (LAMMPS special_bonds 0 0 0), unconfirmed; printed eq. (5) is paper_all_pairs
    seed: int = 12001

    def __post_init__(self):
        if self.n_bonds < 3 or self.r0 <= 0 or self.dt <= 0 or self.gamma <= 0 or self.kbt <= 0:
            raise ValueError("invalid positive physical parameter")
        if self.r0 != 1.0:
            # The printed bending energy kappa/4 sum|b_i-b_(i+1)|^2 (Eq. 4) is written in units of r0; with any other r0
            # the persistence length would be xi_P*r0**3, not xi_P, so the label xi_P/L would be wrong.
            raise ValueError("r0 must be 1: lengths are in units of the bond length r0 (Eq. 4 is written in those units)")
        if self.xi_over_l <= 0 or self.sigma_over_r0 <= 0:
            raise ValueError("xi_over_l and sigma_over_r0 must be positive")
        if self.nonbonded_policy not in POLICIES:
            raise ValueError(f"nonbonded_policy must be one of {tuple(POLICIES)}")
        if self.dt > 0.2 * self.gamma / self.k_s + 1e-18:
            raise ValueError(f"dt={self.dt:g} exceeds conservative bond limit {0.2*self.gamma/self.k_s:g}")

    @property
    def n_beads(self): return self.n_bonds + 1
    @property
    def length(self): return self.n_bonds * self.r0
    @property
    def sigma(self): return self.sigma_over_r0 * self.r0
    @property
    def cutoff(self): return 2.0 ** (1.0 / 6.0) * self.sigma
    @property
    def k_s(self): return self.spring_factor * self.kbt / self.r0**2
    @property
    def kappa(self): return self.xi_over_l * self.kbt * self.length
    @property
    def f_p(self): return self.pe * self.kbt / self.length**2
    @property
    def gamma_l(self): return self.gamma * self.n_beads / self.length
    @property
    def tau(self): return self.length**3 * self.gamma_l / (4.0 * self.kbt)
    @property
    def contour_time(self): return self.length * self.gamma_l / max(self.f_p, 1e-300)
    @property
    def flexure(self): return self.pe / self.xi_over_l
    @property
    def noise_scale(self): return math.sqrt(2.0 * self.kbt * self.dt / self.gamma)
    @property
    def policy_code(self): return POLICIES[self.nonbonded_policy]
    def steps_for_tau(self, units=1.0): return int(math.ceil(units * self.tau / self.dt))
    def canonical(self):
        d = asdict(self); d.update({"length":self.length,"k_s":self.k_s,"kappa":self.kappa,"f_p":self.f_p,"tau":self.tau})
        return d
    def fingerprint(self):
        raw=json.dumps(self.canonical(),sort_keys=True,separators=(",",":"),allow_nan=False).encode()
        return hashlib.sha256(raw).hexdigest()
