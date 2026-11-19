# Visualization note (2026-10-02)

**Not a replication.** Each run below covers 100 simulation time units, which is 3.96e-4 τ (τ = 252 500). That is far below 5 τ, and no transient was discarded. The films and panels are integrity visualizations, not paper figures. Every frame and panel is titled to say so.

Input check: the given zip `spin-chain-paper-v2/spin-chain-paper-v2.zip` has SHA-256 `8137ba53164d3d376971b39edce20f38a5b0dc85e8d330921d295ade04e642fb` (matches). Its contents are identical to the working tree, excluding caches and `output/`.

## What was added
| File | Change |
|---|---|
| `src/spinchain/visualize.py` (new) | Reads a schema-3 HDF5 file and refuses one without `origin_convention` or `positions`. Writes `<name>.gif` (film), `<name>.png` (static panel) and `<name>.json` (summary). |
| `src/spinchain/cli.py` | Added `visualize PATH [--passive P] [--out-dir output/viz] [--fps 12] [--max-frames 200]`. Added `--window-start-tau/--window-end-tau/--window-reason` to `analyze`. `paper-reproduce` and `headless` are unchanged. |
| `tests/test_visualize.py` (new, 4 tests) | Refusal without `origin_convention`; frame cap; GIF frame count equals the number of selected saved frames (no interpolation); the reset check flags a recentered COM. |

Tests: **44 passed**, the 40 existing ones unchanged plus 4 new. The 284 warnings are pyparsing deprecations raised inside matplotlib's mathtext.

**Film.** Lab frame, styled after Fig. 2: body gray, leading tip (code bead N = paper bead 1) black, centre-of-mass trail red. The trail uses every saved frame up to the current one. Each frame prints Pe, ξ_P/L, t/τ and s.
- Frames are real saved configurations, every k-th with k = ceil(saved / 200). There is no interpolation.
- 12 fps was requested. GIF stores delays in 10 ms units, so Pillow writes 80 ms per frame, which plays at 12.5 fps.

**Panel.**
1. s(t)
2. COM path, with 4 saved snapshots
3. |s| histogram (bin 0.25), with the symmetrized β₂ printed as a number only
4. Railway MSD of paper bead j = 25 (code 76) for i = 1, 2, 4, 8, against the passive Pe = 0 bead j. Dotted lines mark i r₀/v.
5. MSR of r_e
6. A text column: COM-reset evidence, policy and dt labels

Lags are in simulation time units. D_r is printed only when `analysis.unit_slope_fit` finds a unit-slope run. Otherwise the panel says "no unit-slope region: D_r not reported".

## Commands (run from the tree root, `PYTHONPATH=src`)
```
python -m spinchain headless --regime polymer       --steps 2000000 --sample-every 20000 --output output/viz/polymer-Pe200.h5 --overwrite
python -m spinchain headless --regime weak_spiral   --steps 2000000 --sample-every 20000 --output output/viz/weak_spiral-Pe1000.h5 --overwrite
python -m spinchain headless --regime strong_spiral --steps 2000000 --sample-every 10000 --output output/viz/strong_spiral-Pe5000.h5 --overwrite
```
The CLI can only make Fig. 2 presets, so the passive references came from the same `run()` with Pe = 0. Each uses N = 100, ξ_P/L = 0.2, seed 12001, 2 000 000 steps, and matches its film's cadence:
```
python -c "from spinchain.paper import make;from spinchain.parameters import Parameters;from spinchain.simulator import run;p=Parameters(**{**make('polymer').__dict__,'pe':0.0});[run(p,2000000,f'output/viz/passive-Pe0-every{se}.h5',sample_every=se,overwrite=True) for se in (20000,10000)]"
```
```
python -m spinchain visualize output/viz/polymer-Pe200.h5       --passive output/viz/passive-Pe0-every20000.h5
python -m spinchain visualize output/viz/weak_spiral-Pe1000.h5  --passive output/viz/passive-Pe0-every20000.h5
python -m spinchain visualize output/viz/strong_spiral-Pe5000.h5 --passive output/viz/passive-Pe0-every10000.h5
```
Each run took about 15 s of CPU (numba, float64).

## Results
| Run | Steps | Saved frames | Film frames (stride) | max \|s\| | COM reset? | COM displacement |
|---|---|---|---|---|---|---|
| polymer, Pe = 200, seed 12001 | 2 000 000 | 101 (every 20 000) | 101 (1) | 0.427 | **no** | 2.58 r₀ |
| weak_spiral, Pe = 1000, seed 12002 | 2 000 000 | 101 (every 20 000) | 101 (1) | 0.725 | **no** | 7.77 r₀ |
| strong_spiral, Pe = 5000, seed 12003 | 2 000 000 | 201 (every 10 000) | 101 (2) | 0.498 | **no** | 44.3 r₀ |
| passive references, Pe = 0 | 2 000 000 each | 101 and 201 | – | – | – | – |

**COM reset check, on every film.** The stored COM equals the mean of the stored lab-frame positions exactly. Zero frames after the start have |COM| < 1e-9; the old recenter wrote about 1e-16 there. The centre of mass was **not reset**.

- **Polymer, Pe = 200.** s stays near 0 (max |s| = 0.43). The railway is **incomplete**: the i = 1 transit r₀/v = 50 falls exactly at the maximum lag of 50, and i = 2, 4, 8 lie beyond it. The panel says so.
  - The D_r unit-slope rule found D = 9.6e-5 on t = 2–5. That is a short-lag crossover, not the diffusive regime: the window is 100, while 1/D_r,a ≈ 1010. Inspect before using; it is not a result.
- **Weak spiral, Pe = 1000.** max |s| = 0.73, so no spiral in this window.
  - Railway transits at 10, 20 and 40 are inside the window, and the j+1, j+2 and j+4 MSDs approach and cross the passive curve near those lags.
  - This is one seed with 101 frames, and curves dipping below the passive line are not significant.
  - No unit-slope MSR region was found, so D_r is not reported.
- **Strong spiral, Pe = 5000.** **No spiral in this window** (max |s| = 0.50 < 1). The film is kept and the run was not extended.
  - The D_r rule found D = 1.7e-4 on t = 9.5–12. This is again the straight-start transient, not the Fig. 8 diffusive regime, so it is not a result.

## Unchanged limits
- `lammps_default` is a hypothesis.
- dt = 5e-5 is unpublished; it is the code guard 0.2 γ/k_S.
- No transient was discarded.
- No morphology classifier is used. "No spiral in this window" is the brief's reporting rule for max |s| < 1, not a regime label.
- One seed per preset.
- No 5 τ run, no three-seed campaign, and no retuning of κ, f_p, dt or noise.

Not a replication.
