# Verdict plan (2026-10-02)

**Not a replication; the 5τ job is still running.**

Tree: `spin-chain-paper-v2/spin-chain-paper-v2 (3).zip`, SHA-256 `617bc22782bb84bcf77c5ee071ecc12a73e5ffd277241caa81a0e0c5b3e9cc84` (matches). Its contents are identical to the working tree, excluding caches and `output/`.
This session wrote the verdict only. Nothing was run on the production job. κ, f_p, dt and noise are unchanged.

## The running job, read-only
Only `pid.txt`, `run.log` and the process table were read. The HDF5 file was not opened, and the job was not stopped, resumed or written to.

| | |
|---|---|
| Job | strong_spiral, paper_all_pairs, seed 12013, 5τ = 25 250 000 000 steps, cadence 0.25 contour |
| `pid.txt` | `launcher PID=9568 interpreter PID=11268 started=2026-10-02T01:38:22.7700703+05:30` |
| Status at 2026-10-02 01:50:14 +05:30 | PID 11268 alive (702 s CPU), launcher 9568 alive |
| Last log line | `progress: step=86,860,000 contours=21.500 tau=0.0172` |
| Observed rate | ≈ 1.2e5 steps/s, giving about 58–62 h in total; expected end ≈ 2026-10-04, 11:00–16:00 +05:30 |
| Directory | `output/reproduction/paper_all_pairs-seed12013/` |

If the process is found dead later: record the last `run.log` line, and resume only in a session that is asked to. The resume command is in CADENCE_NOTE.md.

### Interruption and resume #1 (recorded 2026-10-03)
- **Stop.** The first segment stopped at 2026-10-02 03:38:09 +05:30, with nothing on stderr and no reboot or sleep. 14 s later Windows logged a service install and the "Claude" service being disabled, which suggests a Claude app update killed the process tree. That is inferred from timing, not proven.
  - Last line in `run.log`: `progress: step=862,540,000 contours=213.500 tau=0.1708`.
  - Checkpoint at step 862 540 000.
- **First segment.** Checked read-only after the stop: `strong_spiral-seed12013.h5` is readable, holds 855 frames from step 0 to 862 540 000 at uniform 1 010 000-step cadence, has COM == mean(positions), and is marked `complete=False`. The stitcher accepts that for a non-final segment.
- **Resume.** On 2026-10-03 00:50:29 +05:30 the user started `resume1.cmd` by double-click, in their own session rather than from Claude. It uses the in-tree venv `.venv311` (same pinned versions) because the earlier venv lived in a temporary session folder.
  - Processes: interpreter PID 16532, cmd launcher PID 23184.
  - Output: `strong_spiral-seed12013.from862540000.h5`, logging to `run-resume1.log` and `run-resume1.err`.
  - The physics fingerprint is unchanged. The first resumed progress line was `progress: step=864,560,000 contours=214.000 tau=0.1712`.
- **Remaining.** 24 387 460 000 steps, about 56–68 h at 1.0–1.2e5 steps/s, so the end is expected around 2026-10-05.
- **Completion check.** The job is finished when `run-resume1.log` ends with `"completed": true`, not `run.log`. The stitch glob `strong_spiral-seed12013*.h5` picks up both segments.
- **If it dies again.** Double-click `resume1.cmd` again. It resumes from the newest checkpoint into another `.from<step>.h5` file and appends to the same logs.

## After the job ends
Run from the tree root with `PYTHONPATH=src`, only once PID 11268 has exited and `run.log` ends with `"completed": true`.

1. Stitch the segment or segments (command from CADENCE_NOTE.md):
```
python -c "import glob;from spinchain.stitch import write_stitched;print(write_stitched(glob.glob('output/reproduction/paper_all_pairs-seed12013/strong_spiral-seed12013*.h5'),'output/reproduction/strong_spiral-paper_all_pairs-seed12013-stitched.h5'))"
```
2. Verdict:
```
python tests/verdict.py output/reproduction/strong_spiral-paper_all_pairs-seed12013-stitched.h5 --out output/reproduction/verdict-paper_all_pairs-seed12013.json
```

## What `tests/verdict.py` does
- **Refuses** (exit code 2, no JSON) a file:
  - without `origin_convention`;
  - without `positions`;
  - whose stored COM is not exactly the mean of its stored positions;
  - shorter than 5τ, where duration is the last `time_tau` minus the first;
  - whose saved frames are not uniformly spaced.
- **Window:** drops the first 0.5τ. That is the only selection, and the JSON records it as `"drop the first 0.5 tau; no morphology cut"`. For this job the window is 4.5τ, about 22 500 frames.
- **Reports these and only these**, on the window:
  - `duration_tau`, for the whole run
  - `fraction_abs_s_gt_2`
  - `longest_plateau_contours`: the longest run of consecutive frames with |s| > 2, as t_last − t_first, so a one-frame spike counts 0
  - `beta2_symmetrized`
  - `abs_s_histogram`: bin 0.25
  - `com_msd_loglog_slope`: least-squares slope of log MSD against log lag, from one sample to half the window
  - `msr_unit_slope_region_exists`: the existing `unit_slope_fit` rule, |slope − 1| ≤ 0.1 over ≥ 3 consecutive segments. When it is false, `D_r` is `null`.
- **Provenance** is also recorded: the input path, `stitched_from`, the policy and its hypothesis status, the timestep status, and parameters.
- **Descriptions.** The JSON holds text, not a pass bit, and lists every description whose rule holds; the list may be empty or hold more than one.
  - *predominant*: the |s| > 2 fraction is above 1/2 and β₂ < 2
  - *intermittent*: a |s| > 2 plateau longer than 2 contours exists and the fraction is below 1/2
  - *absent*: no such plateau exists
  - The JSON states that these are descriptions, not a claim that a figure was reproduced.

Tests (`tests/test_verdict.py`, 7):
- one synthetic s(t) series each for predominant, intermittent and absent (single-sample spikes);
- a real 0.008τ integrator file refused as shorter than 5τ;
- a missing `origin_convention` refused;
- a COM edited to differ from mean(positions) refused;
- the window start recorded at 0.5τ on a real 5.5τ file. That file uses a 4-bond chain with τ = 20 time units; it is an integrity size, not physics.

Full suite: **63 passed**. The script also ran end-to-end on a scratch 5.5τ file, writing the JSON, and on a short file, refusing with exit code 2.

## What the verdict will and will not mean
- One seed, one policy, one Fig. 2 point. Even *predominant* would describe this trajectory only. It would not reproduce Fig. 2c, the phase diagram, or the policy choice.
- The gate already showed that seed, not policy, decided formation over 10 contours.
- dt = 5e-5 is the code guard, not the ESI timestep. The special_bonds line is unknown, so lammps_default and paper_all_pairs remain unresolved.

Not a replication; the 5τ job is still running.
