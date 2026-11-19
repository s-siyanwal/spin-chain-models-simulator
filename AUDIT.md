# Scientific and software audit (v2.0)

## Verdict
The repository can test and reproduce the paper's model and qualitative Figure-2 regimes. A quantitative replication is not established until independent trajectories reproduce the reported distributions and curves for at least 5 tau after transient removal. Short laptop demonstrations are not paper replication.

## Corrections made
1. **Nonbonded topology:** the printed sum in Eq. 5 starts at `j > i`, but the authors used LAMMPS. LAMMPS molecular defaults commonly exclude LJ interactions through three bonded neighbors unless `special_bonds` is changed. The primary `lammps_default` policy therefore excludes 1-2, 1-3, and 1-4 WCA pairs. All-pairs and alternative policies remain available for sensitivity testing because the ESI/input script could not be recovered here.
2. **Time scale:** `gamma_l = gamma*(N+1)/L`; hence `tau=L^3*gamma_l/(4*kBT)`. For N=100, r0=gamma=kBT=1, tau=252500 and 5 tau at the conservative default dt=5e-5 is 25,250,000,000 steps per replicate.
3. **Bending force:** implemented exactly from Eq. 4, including the factor `kappa/4`.
4. **Active force:** `fp*b_i` is split equally between bond endpoints and is not normalized.
5. **Observables:** added the paper's bead/COM/orientation foundations and dynamic structure factor. Kurtosis is the standardized fourth moment of the symmetrized spiral-number distribution.
6. **Reproducibility:** counter-based Gaussian noise makes chunking and restart deterministic. Checkpoints are atomic and fingerprinted.
7. **Data integrity:** HDF5 outputs have explicit incomplete/complete state and refuse accidental overwrite.
8. **Evidence:** morphology thresholds from prior versions are retained only as optional regression concepts; they are not substituted for the paper's phase classification.

## Remaining evidence gaps
- The exact timestep and parameter table are stated by the paper to be in the ESI. The public main PDF does not provide each run's timestep.
- The in-house LAMMPS modifications and exact `special_bonds` command were not published in the main article.
- Pixel-perfect figure reproduction requires the original parameter grid, production seeds, sampling cadence, and ESI settings.
- The RX 560X's Vulkan float64 capability and numerical parity must be tested on the target laptop before GPU production. Taichi Vulkan is therefore experimental, not silently selected.

## Hardware decision
N=101 is too small to amortize GPU launch overhead reliably, while billions of serial integration steps dominate. The validated Numba CPU kernel with a skin-based Verlet list is the safe default. Runtime also benchmarks the direct all-pairs reference. RX 560X Vulkan is useful only if an on-device many-replica batched kernel passes force, trajectory-statistics, and long-run stability tests; a single-chain kernel is unlikely to deliver its theoretical throughput.
