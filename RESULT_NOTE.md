# Result note: 5τ trajectory (2026-10-05)

**Predominant on this seed only; not a replication.**

## The run
| | |
|---|---|
| Point | strong_spiral, Pe = 5000, ξ_P/L = 0.2, N = 100 (κ = 20, f_p = 0.5), seed 12013 |
| Policy | paper_all_pairs (printed eq. (5)). This is unresolved against lammps_default: in the gate, both policies gave the same plateau on this seed |
| dt | 5e-5, unpublished; it is the code's stability guard 0.2 γ/k_S, not the ESI timestep |
| Length | 25 250 000 000 steps = **5.0 τ** = 6 250 contours (contour time 202 time units) |
| Completion | `run-resume1.log` ends with `"completed": true` at step 25 250 000 000, 2026-10-05 06:03:49 +05:30. stderr is empty |
| History | Segment 1 ran from step 0 to 862 540 000 and was killed by an apparent Claude app update. Segment 2 resumed from that checkpoint (`resume1.cmd`, started by the user). Both are recorded in VERDICT_PLAN.md |
| Stitch | `output/reproduction/strong_spiral-paper_all_pairs-seed12013-stitched.h5`: 25 001 frames at 0.25-contour cadence, step 0 to 25 250 000 000, `complete=True`. The two segments matched bitwise at the shared frame, so there is one origin |

## Verdict (`output/reproduction/verdict-paper_all_pairs-seed12013.json`)
Window: the first 0.5 τ is dropped, giving 0.5–5.0 τ and 22 501 samples. There is no morphology cut.

| Quantity | Value |
|---|---|
| JSON description(s) | **`["predominant"]`** |
| Fraction of window samples with \|s\| > 2 | **0.978** |
| Symmetrized β₂ | **1.027** |
| Longest \|s\| > 2 plateau | **2007.75 contours** (≈ 1.61 τ, frames 16 969–25 000, t = 3.394–5.0 τ) |
| COM MSD log-log slope (one sample to half the window) | 1.001 |
| Unit-slope MSR region after 0.5 τ | none, so D_r = null |
| Duration | 5.0 τ |

The rule that produced "predominant" is fraction > 1/2 and β₂ < 2. Neither "intermittent" nor "absent" applies, because a plateau longer than 2 contours exists and the fraction is above 1/2. These words describe this trajectory. They are not a claim that a figure was reproduced.

## Figures (all titled "5τ trajectory, not a paper figure")
- `output/reproduction/result/result_panel.png`: s(t) over the full 5 τ with the dropped 0.5 τ shaded and the longest plateau highlighted; the window |s| histogram (from the verdict JSON); the lab-frame COM path.
- `output/reproduction/result/longest_plateau.gif`: 196 real saved frames of the longest plateau only, about every 41st frame of the 8 032, with no interpolation. The box (±24.6 r₀) is re-centred on each frame's COM. Tip black; red COM trail from the plateau start; s and t in contours on each frame.

## What the figures show (descriptive)
- **s(t).** A spiral forms within the first few contours.
  - Long plateaus at |s| ≈ 4: s ≈ +4 until ≈ 1.08 τ, then s ≈ −4 from ≈ 1.24 τ to ≈ 2.4 τ, then s ≈ +4 to the end.
  - Between them, brief break-ups with re-formation in the opposite handedness. Spirals do break up on this seed: a light-red-square behaviour in Fig. 4e terms, not dark-red.
  - The 0.25-contour cadence shows a few one-frame dips (≈ 2.8 τ and ≈ 3.39 τ). They split plateaus under the consecutive-frame rule.
- **The longest plateau is censored.** It is still running at the end of the run, so 2007.75 contours is a lower bound on that spiral's lifetime.
- **|s| histogram.** 97.7% of window samples lie at 3.5 ≤ |s| < 4.75, mode bin 4.0–4.25. 1.95% are at |s| < 1. The shape is a small low-|s| peak and a dominant high-|s| peak, as the text describes Fig. 4c at Pe = 5000. No quantitative comparison was made, because the paper's values are not digitized.
- **COM path.** Dense, diffusive clusters while wound, linked by long directed excursions while elongated. This is qualitatively what the Fig. 2 caption describes as "separated clusters" in the strong-spiral regime. The COM moved 408.5 r₀ net over 5 τ.
- **Dynamics.** COM MSD slope ≈ 1, so translation is diffusive over the window's lags. There is no unit-slope MSR region, so no D_r is reported. The paper also says the long-time MSR ∝ t regime could not be detected for runs with spontaneous break-up.

## What this is not
- **Not a replication.** One seed and one Figure 2 point.
  - The other gated seed, 12003, formed no spiral in 10 contours under either policy, so seed-to-seed behaviour is unknown.
  - Pe = 200 and 1000 were not run at paper length.
  - Figs. 5, 6, 8, 9 and 10 were not computed from this file. The cadence is too coarse for Fig. 6.
  - No phase diagram.
- **Settings still open.** The policy question is open; this run does not discriminate between policies. dt and the special_bonds line remain unknown, since they are not in the main text or the arXiv source.

Predominant on this seed only; not a replication.

## Commands used
Tree root, `PYTHONPATH=src`, interpreter `.venv311`:
```
python -c "import glob;from spinchain.stitch import write_stitched;print(write_stitched(glob.glob('output/reproduction/paper_all_pairs-seed12013/strong_spiral-seed12013*.h5'),'output/reproduction/strong_spiral-paper_all_pairs-seed12013-stitched.h5'))"
python tests/verdict.py output/reproduction/strong_spiral-paper_all_pairs-seed12013-stitched.h5 --out output/reproduction/verdict-paper_all_pairs-seed12013.json
python tests/result_figures.py output/reproduction/strong_spiral-paper_all_pairs-seed12013-stitched.h5 output/reproduction/verdict-paper_all_pairs-seed12013.json
```
Code added this step:
- `visualize.plateau_film` gained `follow` and `trail_from_first` options.
- New script `tests/result_figures.py`.
- One new test for follow mode.

Tests: 64 passed. No integration was started, and κ, f_p, dt and noise are unchanged.
