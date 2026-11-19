$ErrorActionPreference = "Continue"
powercfg /setactive SCHEME_MIN
$env:NUMBA_NUM_THREADS = "8"
$env:OMP_NUM_THREADS = "1"
$env:MKL_NUM_THREADS = "1"
Write-Host "High-performance power plan requested. Keep AC power connected and cooling unobstructed."
Write-Host "In Windows Settings > System > Display > Graphics, add .venv\Scripts\python.exe and choose High performance (RX 560X)."
