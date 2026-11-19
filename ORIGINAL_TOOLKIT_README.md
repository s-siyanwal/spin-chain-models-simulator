# Spin-chain paper reproduction toolkit v2.0

Audited implementation of R. E. Isele-Holder, J. Elgeti and G. Gompper, *Soft Matter* **11**, 7181–7190 (2015), DOI 10.1039/c5sm01683e.

## What this can claim
It implements the published 2D overdamped Brownian model, generates qualitative polymer/weak-spiral/strong-spiral tests, and supports a checkpointed multi-seed reproduction campaign. It does **not** claim that a smoke test or one short trajectory reproduces the paper.

## Install on DESKTOP-C7IVGS2
Open PowerShell in the extracted folder:

```powershell
powershell -ExecutionPolicy Bypass -File .\setup.ps1
.\.venv\Scripts\Activate.ps1
python main.py hardware
python main.py paper-audit --regime polymer --tau 5 --replicates 3
```

Use AC power, Windows Best performance, and an elevated cooling profile. The Ryzen 5 2500U has four cores/eight threads; the production engine uses a compiled Numba loop and avoids oversubscribing BLAS. Run `tools\windows-performance.ps1` if desired.

## Reproduction
The primary hypothesis is LAMMPS default topology exclusion:

```powershell
python main.py paper-reproduce --regime polymer --tau 5 --replicates 3 --policy lammps_default --i-know
```

Repeat for `weak_spiral` and `strong_spiral`. Also run the all-pairs sensitivity:

```powershell
python main.py paper-reproduce --regime strong_spiral --tau 5 --replicates 3 --policy paper_all_pairs --i-know
```

These jobs are extremely long. `paper-audit` benchmarks the laptop and reports an estimate before launch. Resume support is available to the lower-level `headless` command.

## Short validation

```powershell
python main.py headless --regime polymer --steps 100000 --sample-every 5000 --output output\polymer-smoke.h5
python main.py analyze output\polymer-smoke.h5 --structure
```

## GPU policy
The optional Taichi package can probe Vulkan:

```powershell
pip install -r requirements-gpu.txt
python main.py vulkan-probe
```

Vulkan is not approved automatically. Taichi's Vulkan float64 support is extension-dependent, and the RX 560X must pass a local parity benchmark. Force `.venv\Scripts\python.exe` to **High performance** in Windows Graphics settings to prefer the discrete GPU. Never infer the selected adapter from marketing names alone.

## Physics
- `N+1` beads, `N` bonds, `L=N*r0`
- harmonic bonds with `kS=4000 kBT/r0^2`
- vector bending energy `kappa/4 sum |b_i-b_(i+1)|^2`
- WCA exclusion, free boundaries, no hydrodynamics
- tangential bond force density split between endpoints
- Euler-Maruyama with explicit `gamma` and thermal displacement
- `xiP/L=kappa/(kBT L)`, `Pe=fp L^2/kBT`

See `AUDIT.md` for unresolved ESI and LAMMPS-topology issues.

## Hardware-specific optimization
The runtime benchmarks both direct all-pairs and skin-based Verlet engines, then reports the faster one. The headless simulator uses the Verlet path, which avoids checking every WCA pair each step. On the build host this was roughly 3.1 times faster; treat that only as an implementation check and run `paper-audit` on DESKTOP-C7IVGS2 for the real estimate. The 32 GB RAM is ample; a 10,001-frame position trajectory with 101 beads is about 16 MB before HDF5 overhead/compression. Keep sample cadence moderate to protect the 238 GB SSD.

For independent replicas, run separate processes only after benchmarking one versus two concurrent jobs. The 2500U is a 15 W mobile four-core CPU, so four simultaneous jobs will usually throttle; one or two long jobs on AC power are the sensible starting points.
