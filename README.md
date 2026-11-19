# spin-chain-models-simulator

Brownian-dynamics simulation of a **self-propelled worm-like filament** (an active polymer chain), written as an audited reproduction toolkit for

> R. E. Isele-Holder, J. Elgeti, G. Gompper, *Self-propelled worm-like filaments: spontaneous spiral formation, structure, and dynamics*, Soft Matter **11**, 7181–7190 (2015), DOI [10.1039/c5sm01683e](https://doi.org/10.1039/c5sm01683e).

**Status: this is a model implementation and an in-progress reproduction attempt. It is not a completed replication of the paper.** See [What has and has not been shown](#what-has-and-has-not-been-shown).

> Provenance note: the code in this repository is the "v2.0" toolkit whose notes are dated 2026-09-30 to 2026-10-05. It is not the original 2020 summer-project code; that code is not in this repository.

## What the code does

Python package `spinchain` (`src/spinchain/`), 2D overdamped Brownian dynamics, no hydrodynamics:

- Chain of `N+1` beads with `N` bonds, `L = N·r0`; harmonic bonds (`kS = 4000 kBT/r0²`), vector bending energy `κ/4 Σ|b_i − b_(i+1)|²`, WCA excluded-volume repulsion, free boundaries.
- Active force: a tangential force density along each bond, split between its two end beads, of strength `f_p` (the paper's Eq. 7). Control parameters are `Pe = f_p L²/kBT` and `ξ_P/L = κ/(kBT L)`.
- Euler–Maruyama integration with a compiled **Numba** kernel, a skin-based Verlet neighbour list, and counter-based Gaussian noise so that chunked and resumed runs are bitwise deterministic. Checkpoints are atomic and fingerprinted; HDF5 output (schema 3) stores positions in the lab frame with no recentring, plus an unwrapped end-to-end angle.
- Three presets from the paper's Figure 2 (`polymer` Pe=200, `weak_spiral` Pe=1000, `strong_spiral` Pe=5000, all ξ_P/L=0.2, N=100).
- Four selectable non-bonded exclusion policies (`paper_all_pairs`, `exclude_12`, `exclude_12_13`, `lammps_default`). Which one the authors used is unresolved (see [`AUDIT.md`](AUDIT.md)); the default `lammps_default` is labelled a hypothesis.
- Analysis (reads only the HDF5 file): spiral number and its symmetrized kurtosis β₂, end-to-end distance vs Kratky–Porod, centre-of-mass MSD, mean-square rotation, relative-bead MSD, and the dynamic structure factor.
- Optional Taichi/Vulkan backend (`requirements-gpu.txt`), experimental and not selected automatically.
- Command line: `python main.py {benchmark,paper-audit,headless,paper-reproduce,analyze,visualize,hardware,vulkan-probe}`.

The time step `dt = 5e-5` is the code's stability guard (`0.2·γ/k_S`). It is **not** the paper's value, which is given only in the paper's ESI.

## Install and test

Python 3.11 (the pins in `requirements.txt` target it).

```bash
pip install -r requirements.txt
PYTHONPATH=src python -m pytest -q          # 64 tests
```

Short integrity run (a few thousand steps, about a second; the output goes to the git-ignored `output/`):

```bash
PYTHONPATH=src python tests/smoke_integrity.py --regime polymer --steps 4000 --sample-every 100
```

For a longer headless run and its analysis:

```bash
PYTHONPATH=src python main.py headless --regime polymer --steps 100000 --sample-every 5000 --output output/polymer-smoke.h5
PYTHONPATH=src python main.py analyze output/polymer-smoke.h5 --structure
```

A full 5 τ run is 25 250 000 000 steps (`paper-audit` prints an estimate; the 5 τ run described below took several days on a laptop). [`ORIGINAL_TOOLKIT_README.md`](ORIGINAL_TOOLKIT_README.md) is the toolkit's own README with Windows-specific setup (`setup.ps1`).

## What has and has not been shown

From the toolkit's own notes (read them for details): [`REPRO_NOTE.md`](REPRO_NOTE.md), [`CADENCE_NOTE.md`](CADENCE_NOTE.md), [`PRODUCTION_NOTE.md`](PRODUCTION_NOTE.md), [`VERDICT_PLAN.md`](VERDICT_PLAN.md), [`RESULT_NOTE.md`](RESULT_NOTE.md), [`AUDIT.md`](AUDIT.md).

- One strong-spiral (Pe=5000) trajectory, seed 12013, policy `paper_all_pairs`, was run for 5.0 τ (25 250 000 000 steps). Over the 0.5–5.0 τ window the fraction of samples with |s| > 2 is 0.978, β₂ = 1.027, and the longest plateau is 2007.75 contours (still running at the end, so a lower bound). The notes describe this as "predominant on this seed only".
- That is **one seed at one point**. Another gated seed (12003, `lammps_default`) formed no spiral in 10 contours, so seed-to-seed behaviour is unknown. Pe=200 and 1000 were not run at paper length, and no phase diagram or other paper figure was reproduced quantitatively.
- Policy, time step and the `special_bonds` setting remain open; the run does not discriminate between policies.
- The remaining short runs are integrity tests, not results.

Small result files from these runs are kept in `results/` (verdict JSON, gate JSON/logs/plots, result panel). The notes refer to the original `output/` layout; the large trajectories are not in git (see below).

## Repository layout

| Path | Contents |
|---|---|
| `src/spinchain/` | simulator, forces, engine, RNG, analysis, observables, I/O, CLI |
| `tests/` | 64 pytest tests plus the gate, verdict, figure and smoke scripts |
| `config/benchmarks.yaml`, `main.py`, `setup.ps1`, `tools/` | configuration, entry point, Windows helpers |
| `docs/REPLICATION_PLAN.md`, `*_NOTE.md`, `AUDIT.md` | plans and audit notes |
| `reports/` | LaTeX sources and PDFs of a "beginner" and an "expert" write-up, with figures |
| `results/` | small outputs from the gate and 5 τ runs |

## Not included (left out of git)

- HDF5 trajectories: the stitched 5 τ file (about 33 MB), its raw segment (about 32 MB), smoke and gate trajectories, `.npz` checkpoints and `.gif` animations (largest about 7 MB). They stay with the original working tree.
- The original paper PDF and its arXiv source (copyrighted; fetch from the DOI above).
- An earlier prototype (`active-polymer-simulator` v1.1–v1.5) and AI-prompt/draft files from the same folder.
- The toolkit's original `LICENSE`, `CITATION.cff` and `MANIFEST.json`. See [Authorship](#authorship).

## Authorship

This toolkit is Shivanshu Siyanwal's own work, written in 2026. The toolkit's original `LICENSE` and `CITATION.cff` files, which are not in this repository, named Akanksha Singh; that name is not an authorship claim. The author is credited in [`CITATION.cff`](CITATION.cff), `pyproject.toml`, [`LICENSE`](LICENSE) and on the title pages of the two reports in `reports/`.

## License

MIT, see [`LICENSE`](LICENSE).
