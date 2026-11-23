# Replication plan

## Tier 1: equations and invariants
Run `python -m pytest -q`. This validates force gradients, net active force, topology policies, counter-RNG continuation, HDF5 completion, spiral number, kurtosis, Kratky-Porod, MSD/MSR primitives, and exact duration accounting. Then run `python main.py wlc-check` (about 100 s), which checks the passive limit (Pe = 0) against the exact bending-angle distribution and tangent correlation of the discrete worm-like chain.

## Tier 2: qualitative Figure 2
For Pe 200, 1000, 5000 at xiP/L=0.2 and N=100, run at least three independent seeds. Verify respectively: low winding; intermittent |s| around 2-3; predominant long-lived |s| around 4. These are qualitative checks, not the phase diagram.

## Tier 3: structural result
At Pe=0 and non-spiral active cases, compare mean squared end-to-end distance with Eq. 15. Quantify uncertainty by independent seeds/block means.

## Tier 4: dynamics
Reproduce: railway bead MSD (Fig. 6); diffusive elongated-state MSR and ballistic spiral-state MSR (Fig. 8); COM-MSD crossover and spiral suppression (Fig. 9); dynamic structure factor near qL=5*pi (Fig. 10).

## Tier 5: phase map
Recover the ESI's complete grid and exact timestep table. Run >5 tau after transient removal. Use the symmetrized standardized kurtosis and identify non-breaking spirals exactly as described by the paper. Report confidence intervals and both nonbonded-policy sensitivities until the original `special_bonds` line is known.
