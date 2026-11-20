# Changes after the dated notes

The notes (`*_NOTE.md`, `AUDIT.md`, `VERDICT_PLAN.md`) record the toolkit as of 5 October 2026 and are kept as written. This file lists what changed afterwards.

## November 2026: portfolio pass

### Fixed

- `paper.audit()` (the `paper-audit` command) reported `paper_faithful: True` for `lammps_default`, which `AUDIT.md` and `REPRO_NOTE.md` call a hypothesis. It now reports the policy and its status string. `REPRO_NOTE.md` still describes the old behaviour.
- `Parameters` refuses `r0 != 1`. The printed bending energy is written in units of r0; with r0 = 2 a simulated chain was r0³ = 8 times stiffer than ξ_P/L says (measured ℓ_p/ξ_P = 7.99).
- `simulator.run()` checks its cadence before opening any file. Before the fix, `chunk=0` hung, `sample_every=0` crashed after creating the output, and a negative `sample_every` saved frames off the recorded cadence.
- `CADENCE_NOTE.md` and `PRODUCTION_NOTE.md` embedded images from the git-ignored `output/`. They now point to the byte-identical copies in `results/gate/`.

The force kernel, integrator, noise, time step and default policy are unchanged. Physics fingerprints do not change, so existing checkpoints and HDF5 files remain valid.

### Added

- A passive-limit check of the bending term against the exact discrete worm-like-chain statistics: `src/spinchain/validation.py`, `python main.py wlc-check`, and the figure and data in `results/validation/`. The new observables `turning_angles` and `tangent_correlation` come with analytic tests.
- `notebooks/reproduce.ipynb` and `requirements-notebook.txt`. A second CI job executes the notebook.
- More tests: 64 → 127, covering the force law, the noise and integrator, determinism and integrity refusals, the observables, and the verdict rules.
- `CITATION.cff`. The author is now credited in `pyproject.toml`, `LICENSE`, the README and the reports' title pages; the report PDFs were rebuilt with their content unchanged.
- A rewritten README.
