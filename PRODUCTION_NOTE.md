# Production note (2026-10-02)

> **Superseded by [CADENCE_NOTE.md](CADENCE_NOTE.md).** The seed-12003, 0.01 τ job described below was stopped at step 78 780 000 and will not be resumed. Its folder is now `output/reproduction/STOPPED-paper_all_pairs-seed12003-wrong-seed/`. The unconfounded gate showed that formation follows the seed, not the policy.

**Started, not finished; not a replication.**

Tree: `spin-chain-paper-v2/spin-chain-paper-v2.zip`, SHA-256 `24768cadcd3ed581336e26890d4066fad7e6d85b5a268e640cd7ada949c4c925` (matches). Its contents are identical to the working tree, excluding caches and `output/`.

## Formation gate
Setup: strong_spiral (Pe = 5000, ξ_P/L = 0.2, N = 100, κ = 20, f_p = 0.5), dt = 5e-5 (unpublished), schema 3, lab frame, recenter_every = 0.
- One contour time is L·γ_l/f_p = 100 × 1.01 / 0.5 = 202 time units, or 4 040 000 steps.
- 10 contours is **40 400 000 steps** (0.008 τ).
- `sample_every` was 101 000 steps (0.025 contour), giving 401 saved frames with positions.
- Checkpoint and angle sidecar were written at **every** saved frame. That is stricter than the every-50-frames requirement, and checkpoint cost is negligible at this cadence.

**Pass rule, the only test:** a run of consecutive saved frames with |s| > 2 lasting at least 2 contour lengths (t_last − t_first, so a one-frame spike scores 0), **and** the stored COM equals the mean of the stored positions with no reset. The rule is implemented in `observables.longest_run_above`, which has a unit test, and evaluated in `tests/formation_gate.py`. No other morphology measure is used.

![gate s(t)](output/gate/gate_s_t.png)

| Policy | Seed | Steps | Contours | max \|s\| | Longest \|s\|>2 run | COM == mean(positions) | Frames with \|COM\|<1e-9 after start | COM displacement | Gate |
|---|---|---|---|---|---|---|---|---|---|
| lammps_default | 12003 | 40 400 000 | 10.0 | 1.68 | none | yes | 0 | 197.4 r₀ | **failed** |
| paper_all_pairs | 12013 | 40 400 000 | 10.0 | 4.44 | **3.80 contours** (153 frames, t = 4.8 → 8.6 contours, s ≈ −4.1) | yes | 0 | 258.8 r₀ | **passed** |

**Passing policy: `paper_all_pairs`** (printed eq. (5), all pairs j > i).
- The plateau is clean in s(t): it forms at 4.8 contours, holds at s ≈ −4.1 for 3.8 contours, then breaks up.
- The COM trail in the film flattens while the chain is wound, consistent with the paper's "almost purely diffusive" spiral state.

Caveats:
- **Policy and seed are confounded.** The gate compared lammps_default/12003 with paper_all_pairs/12013, as specified. One spiral in one seed does not show that lammps_default cannot spiral, or that paper_all_pairs always does.
- **The production seed was not itself gated.** Production uses paper_all_pairs with seed **12003**, a seed that was gated only under lammps_default.
- **lammps_default is still a hypothesis.** paper_all_pairs passing is not evidence that it is the paper's policy; the original special_bonds line is still unavailable.

Gate commands (run from the tree root with `PYTHONPATH=src`). Each took about 5–6 min, run in parallel:
```
python -m spinchain headless --regime strong_spiral --policy lammps_default --steps 40400000 --sample-every 101000 --output output/gate/strong_spiral-lammps_default-seed12003.h5 --checkpoint output/gate/strong_spiral-lammps_default-seed12003.ckpt.npz --overwrite
python -c "from spinchain.paper import make;from spinchain.parameters import Parameters;from spinchain.simulator import run;p=Parameters(**{**make('strong_spiral','paper_all_pairs').__dict__,'seed':12013});b='output/gate/strong_spiral-paper_all_pairs-seed12013';print(run(p,40400000,b+'.h5',sample_every=101000,checkpoint=b+'.ckpt.npz',overwrite=True))"
python tests/formation_gate.py output/gate/strong_spiral-lammps_default-seed12003.h5 output/gate/strong_spiral-paper_all_pairs-seed12013.h5
```
The CLI has no `--seed` flag, so the seed-12013 run used the same `run()` from Python.

Outputs in `output/gate/`:
- `gate_s_t.png` and `gate_result.json`
- for each run, a film `.gif` (134 real saved frames, stride 3, 80 ms per frame) and a panel `.png`. All are titled "formation gate, not a paper figure", in the lab frame, with tip black and COM trail red.
- The films keep a fixed lab-frame box around the whole trajectory (COM travels about 250 r₀), so the wound spiral looks small.

## Code changes this session
- `cli.py`, `paper-reproduce` only:
  - `--sample-every` (default **0.01 τ** = 50 500 000 steps at dt = 5e-5, replacing `steps//10000`).
  - `--resume` continues from the checkpoint into a **new segment file** `<regime>-seed<seed>.from<step>.h5`. It runs only the remaining steps and never overwrites earlier frames.
  - Prints a plan line: total, done and remaining steps, cadence, output, checkpoint.
  - Helper `reproduce_sample_every`.
- `simulator.py`: one line, `log.f.flush()` before each checkpoint, so frames named by a checkpoint are on disk.
- `observables.py`: `longest_run_above` (the plateau rule).
- `visualize.py`: `title` parameter.
- New tests:
  - `tests/test_reproduce_cli.py` (3): default 0.01 τ; resume is bitwise equal to an uninterrupted run; the earlier segment is kept.
  - One `longest_run_above` test in `tests/test_figures.py`.
- New script `tests/formation_gate.py`.
- Tests: **48 passed**.

## The 5τ job
| | |
|---|---|
| Regime / policy / seed | strong_spiral / **paper_all_pairs** / **12003** |
| Steps | **25 250 000 000** (5 τ, τ = 252 500; 1 250 contours) |
| Sample interval | 50 500 000 steps = 0.01 τ, so 501 frames with positions; checkpoint at every frame |
| Benchmark (20 000 steps, `spinchain benchmark --regime strong_spiral --policy paper_all_pairs`) | 115 335 steps/s (Verlet), giving **≈ 60.8 h (2.5 days)** wall clock |
| Rate observed in the gate (with spirals) | 120 380 steps/s, giving ≈ 58 h |
| Started | 2026-10-02 01:19:52 +05:30, detached. Launcher PID 12856, interpreter PID 21328 |
| Output | `output/reproduction/paper_all_pairs/strong_spiral-seed12003.h5` |
| Checkpoint | `output/reproduction/paper_all_pairs/strong_spiral-seed12003.h5.checkpoint.npz` (+ `.angle.json`) |
| Logs | `output/reproduction/paper_all_pairs/run.log`, `run.err` |

Launch command (from the tree root, `PYTHONPATH=src`):
```
python -m spinchain paper-reproduce --regime strong_spiral --policy paper_all_pairs --tau 5 --replicates 1 --output-dir output/reproduction/paper_all_pairs --i-know
```
**Resume command**, after an interruption (from the tree root, `PYTHONPATH=src`):
```
python -m spinchain paper-reproduce --regime strong_spiral --policy paper_all_pairs --tau 5 --replicates 1 --output-dir output/reproduction/paper_all_pairs --i-know --resume
```
This continues from the last checkpoint into `strong_spiral-seed12003.from<step>.h5`. Keep both segment files. Their frames meet at the checkpoint step, which appears in both files, and the end-to-end-angle branch continues through the sidecar.
- If the job dies before its first checkpoint (≈ 7 min), there is nothing to resume. Rerun the launch command with `--overwrite`.
- To stop the job: `Stop-Process -Id 21328`.
- Do not open the running `.h5` with h5py while the job writes it; HDF5 file locking can interfere.

## Still open
- **Sampling.** The 0.01 τ cadence is 12.5 contours, coarser than the 3.8-contour plateau seen in the gate.
  - The 501 frames sample p(|s|), β₂, ⟨r_e²⟩ and long-lag COM MSD and MSR.
  - They cannot resolve individual plateaus, Fig. 6 railway lags or the short-lag part of Fig. 10.
  - A denser `--sample-every` can be chosen at the next resume. It applies only to that new segment.
- **Segment stitching.** `analyze()` reads one file; stitching segments is not implemented.
- **Unpublished inputs.** dt = 5e-5 is unpublished. The ESI timestep and the special_bonds line are unknown.
- **No transient discard yet.** Analysis must set `window_tau` explicitly.
- **Scope.** No other regime, seed or policy was launched.

Started, not finished; not a replication.
