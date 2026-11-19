# Repro note: origin fix and Figure 2–10 inputs (2026-10-02)

**Nothing here is a replication.** Every run below is far shorter than 5 τ after a transient. These are integrity tests only. No figure has been reproduced.

Archive verified before editing: `spin-chain-paper-v2.zip` SHA-256 `9ecd7cbc97973bed717b182185f6b88b9409848d24ee52a679a298635d9a0514` (matches).
Environment: Python 3.11.15, numpy 2.1.3, numba 0.61.2, h5py 3.13.0, pytest 8.3.5 (pinned `requirements.txt`).
Sources read: `c5sm01683e.pdf` (main text), plus the arXiv LaTeX `arXiv-1508.01418v1.tar.gz` (`manuscript.tex`), which the user supplied. The LaTeX matches the PDF on every equation used here, and it contains no timestep or `special_bonds` information.

## Files changed
| File | Change |
|---|---|
| `src/spinchain/engine.py` | `recenter_every` defaults to 0 in `advance` and `advance_verlet`. Recentering is opt-in, needs an `offset` accumulator, and raises without one. The subtracted COM is added to `offset`. RNG, Verlet rebuild rule and force kernels are untouched. |
| `src/spinchain/simulator.py` | `run()` (used by `headless` and `paper-reproduce`) calls `advance_verlet(..., recenter_every=0)`. It unwraps the end-to-end angle on a global grid every `ANGLE_EVERY=1000` steps, independent of chunk size and resume point, and carries the angle branch across resume in `<checkpoint>.angle.json`. It also records the largest unwrap increment. |
| `src/spinchain/io.py` | HDF5 schema 3. New datasets: `time` (simulation units) and `ee_angle` (unwrapped). New attrs: `origin_convention`, `leading_tip`, `nonbonded_policy`, `nonbonded_policy_status`, `timestep_status`, `sample_every`, `start_step`, `ee_angle_branch`, `ee_angle_unwrap_every_steps`, `ee_angle_max_increment_rad`. |
| `src/spinchain/observables.py` | Added `unwrap_next`, `relative_bead_msd` (Fig. 6), `unit_slope_fit` (Fig. 8 D_r rule), `structure_decay_time` (eq. 29), `abs_spiral_histogram` (Fig. 4c), `active_msd_eq25` (comparison curve). Existing kernels are unchanged. |
| `src/spinchain/analysis.py` | Rewritten as the "second program", which reads only the HDF5 file. Details below. |
| `src/spinchain/parameters.py` | One line: default policy `exclude_12_13` → `lammps_default`, labelled hypothesis. |
| `tests/test_origin.py` (new, 8 tests) | Origin regression, chunk and resume determinism, offset accumulator. |
| `tests/test_figures.py` (new, 13 tests) | Analytic and synthetic checks for every new observable, the leading-tip map, window, schema refusal, and policy-default agreement. |
| `tests/smoke_integrity.py` (new, not collected by pytest) | Smoke driver. |

Not changed: `cli.py`, `paper.py`, `forces.py`, `rng.py`, `checkpoint.py`. `cli.py` needed no edit because the CLI and `make()` already defaulted to `lammps_default`.

Tests: **40 passed** (19 original + 21 new). The origin regression `test_stored_com_not_reset_on_recenter_boundary` was run against the untouched v2 tree first and **failed** there: stored COM norm was 2.0e-17 at step 10000. It passes now.

## Origin convention
- The frame is the lab frame. The origin is the centre of the initial straight chain. Positions are never recentered, and there are no periodic images.
- `com` is the mean of the stored lab-frame positions. It is unwrapped, and a test asserts exact equality with the positions.
- Every frame in a file shares one origin. A test checks the stored positions bitwise against an independent un-recentered integration.
- A resumed file continues the same origin, because the checkpoint holds lab-frame positions. A test checks resume bitwise against an uninterrupted run, including `ee_angle`.
- Float64 is enough for this: at expected excursions of about 10³ r0, rounding is about 1e-13 r0.
- **Old schema-2 files cannot be repaired.** They store recentered frames and no offset. `analyze()` refuses them.

## Index convention (relevant to Figure 6)
Eq. (7) is `F_p = Σ f_p r_{i,i+1}` with `r_{i,j} = r_i − r_j`, so in the paper bead 1 leads. The Fig. 6 caption confirms this: "leading bead j". The code kernel applies `+f_p (r_{i+1} − r_i)`, so code bead `n_beads−1` leads. The two are the same physics under the relabelling i → N−i, so the force was **not** changed. The mapping is paper bead k ↔ code index `n_beads − k` (`analysis.paper_bead_to_code`, tested against the sign of the net propulsion force). For N = 100, paper bead j = 25 is code bead 76, and bead j+i is code bead 76−i. The sign of s and of r_e flips under this relabelling. |s|, β₂, the MSR, the COM MSD and S(q,t) are unaffected.

## What `analysis.analyze()` produces
- **Window:** `window_tau=(start, end)` is explicit, and `window_reason` is recorded. The default is `(0, None)` with a warning that the transient is included. Nothing discards samples automatically. For a run that never leaves a spiral, a human sets the window to the spiral state after inspecting s(t). No regime is inferred from β₂, |s| or r_e.
- **Fig. 4:** symmetrized β₂ (unchanged kernel), median |s|, and the |s| histogram (bin 0.25 turn, edges, counts, density).
- **Fig. 5:** ⟨r_e²⟩/L² and Kratky–Porod eq. (15).
- **Fig. 9:** COM MSD with lags in time units and in τ.
- **Fig. 8:** MSR from the stored unwrapped angle, not re-unwrapped at sample cadence. D_r is fitted only on the longest contiguous run of ≥3 log-log segments with |slope−1| ≤ 0.1. Otherwise the result is `None`.
- **Comparison values, not fit targets:** v_c = f_p/γ_l, D_t = kT/(γ_l L), D_r,a = v_c/ξ_p (eq. 23), the rod bound D_r,p = 9/(4τ), and v₀ = f_p√⟨r_e²⟩/(γ_l L) (eq. 27). Eq. (24) is flagged as valid for the polymer regime only.
- **Fig. 6:** MSD of bead j+i relative to bead j, for i ∈ {1, 2, 4, 8}. The paper does not list its i values. Also the self-MSD of bead j, and the passive reference from a Pe = 0 file (`passive_path=`, which must have the same N and sample interval). There is no pass/fail classifier.
- **Fig. 10:** S(q,t)/S(q,0) at q = 5π/L, averaged over 16 directions, using the existing eq. (28) kernel on single-origin positions. τ_S comes from eq. (29), and **its limits are mine**: t = 0 to t_max = ⌊window samples / 2⌋ sample intervals, trapezoid rule on the reported lag grid. The paper writes a bare ∫_t. The value S/S0 at t_max is reported so truncation is visible.
- **Lags** are always reported in time units (`lag_time`), never as sample indices.

## Unpublished timestep
dt = 5e-5 is the code's stability guard 0.2 γ/k_S. It is **not** a published value: the paper says only that Δt "was adjusted… to ensure stable simulations", with the list in the ESI. Every trajectory records `timestep_status` saying so. Placeholder: `ESI_TIMESTEP_PER_RUN` (unknown).

## Two policies still to be compared
- `lammps_default` (excludes separations 1, 2, 3) is now the default of `Parameters`, `make()` and the CLI. A test asserts that all three agree. It is a **hypothesis**: the LAMMPS `special_bonds lj 0 0 0` default is assumed, and the original input is unavailable. Placeholder: `ESI_SPECIAL_BONDS_LINE` (unknown).
- `paper_all_pairs` (printed eq. (5), all j > i) remains one flag away: `--policy paper_all_pairs`.
- `paper.audit()` still prints `paper_faithful = (policy == 'lammps_default')`. That is not evidence; `paper.py` was outside the edit scope.

## Smoke results (not a replication)
Exact command, the few-thousand-step smoke requested (Pe = 200, ξ_p/L = 0.2, N = 100, seed 12001, plus a Pe = 0 reference; about 1 s per run):
```
python tests/smoke_integrity.py --regime polymer --steps 4000 --sample-every 100
```
- Unwrapped COM goes from (0, 6.3e-9) to (−0.0160, 0.1259), |Δ| = 0.127, in the lab frame with no recentering. Over 0.2 time units this is thermal, not ballistic: v₀t ≈ 0.004 against √(4D_t t) ≈ 0.09.
- Spiral number: max |s| = 0.137, median |s| = 0.054, so near 0 as expected at Pe = 200 and short time.
- Railway: **the touch cannot be seen in 4000 steps.** Bead speed is f_p r₀/γ = 0.02 r₀ per time unit, so i = 1 needs about 50 time units (~10⁶ steps), and the max lag here is 0.1. MSD₍ⱼ₊ᵢ₎/MSD_passive stays ≥ 13.

Supplementary integrity run, needed to show the railway touch (weak_spiral preset Pe = 1000, ξ_p/L = 0.2, seed 12002, 100 time units ≈ 4e-4 τ, about 18 s per run):
```
python tests/smoke_integrity.py --regime weak_spiral --steps 2000000 --sample-every 1000
```
- The minimum of MSD₍ⱼ₊ᵢ₎/MSD_passive falls at lag 14.8, 21.9 and 40.2 for i = 1, 2, 4 (predicted transit i r₀/v = 10, 20, 40), with ratios 0.98, 0.88 and 0.69. For i = 8, transit is not reached within the window.
- This is one seed with few independent time origins at long lag, so the ratios below 1 are not significant. Fig. 6 itself uses ξ_p/L = 0.3, not 0.2.
- COM |Δ| = 7.77. Max |s| = 0.85. The largest unwrap increment of the end-to-end angle is 0.0195 rad, far below π.
- The D_r unit-slope rule returned D = 1.27e-4 on lags 1.15–1.75. The window (100) is shorter than 1/D_r,a ≈ 200, so this is a short-lag crossover passing through slope 1, **not** the diffusive regime. The rule is mechanical, and a human must look at the MSR curve before using a fitted D_r.
- S/S0 = 0.94 at the upper limit, so τ_S here (24.8) is truncated and meaningless.

CLI path (also checked):
```
python -m spinchain headless --regime polymer --steps 4000 --sample-every 100 --output output/smoke/cli-polymer-4000.h5 --overwrite
python -m spinchain analyze output/smoke/cli-polymer-4000.h5 --structure --json-out output/smoke/cli-polymer-4000.json
```
(Run from the repo root with `PYTHONPATH=src`, or after `pip install -e .`.) The CLI `analyze` cannot pass `window_tau` or `passive_path`, because `cli.py` was outside the edit scope. Use the Python API for those.

## Paper vs. brief, recorded rather than reopened
- In both `c5sm01683e.pdf` and `manuscript.tex`, eq. (8) prints ξ_p/L = κ/(k_BT L), eq. (3) prints (k_S/2)(|r|−r₀)² with no square root, eq. (1) has γṙ on the left, and the text gives D_t = k_BT/(γ_l L). These are the forms the code already uses. The "printed" variants described in the brief do not appear in these copies, so there is no code consequence.
- The Fig. 3 caption in both sources prints (a,b) ξ_p/L = 0.4, Pe = 10 000, 𝔉 = 25 000 and (c) ξ_p/L = 1.0, Pe = 900, 𝔉 = 900. Both agree with 𝔉 = Pe/(ξ_p/L). The 𝔉 = 2500 label described in the brief is not in these copies. The code follows eq. (10) either way; nothing was changed.

## Remaining limits on quantitative replication
- The ESI timestep (`ESI_TIMESTEP_PER_RUN`) and the original special_bonds line (`ESI_SPECIAL_BONDS_LINE`) are both missing.
- The length of the "initial period" the paper discards is unspecified.
- The integration limits of eq. (29) are unspecified.
- The i values in Fig. 6 are unspecified.
- `paper-reproduce` still samples every `steps//10000` steps (about 126 time units at 5 τ), set in `cli.py`, which was outside scope. The frames are no longer recentered, but that cadence is too coarse for Fig. 6 and Fig. 10 lags and for the τ_S integral. Choose `sample_every` with `headless` before any production run.
- 5 τ is 2.525e10 steps per replicate at dt = 5e-5. It was not launched.
