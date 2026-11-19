$ErrorActionPreference = "Stop"
$py = Get-Command py -ErrorAction SilentlyContinue
if (-not $py) { throw "Python launcher not found. Install 64-bit Python 3.11." }
$ver = & py -3.11 -c "import sys; print(f'{sys.version_info.major}.{sys.version_info.minor}')"
if ($ver -ne "3.11") { throw "Python 3.11 x64 is required." }
if (-not (Test-Path .venv)) { & py -3.11 -m venv .venv }
& .\.venv\Scripts\python.exe -m pip install --upgrade pip setuptools wheel
& .\.venv\Scripts\python.exe -m pip install -r requirements.txt
& .\.venv\Scripts\python.exe -m pip install -e .
New-Item -ItemType Directory -Force logs,output | Out-Null
& .\.venv\Scripts\python.exe main.py hardware | Out-File -Encoding utf8 logs\hardware.json
& .\.venv\Scripts\python.exe -m pytest -q
& .\.venv\Scripts\python.exe main.py benchmark --steps 20000 | Out-File -Encoding utf8 logs\benchmark.json
Write-Host "CPU production environment ready. Optional experimental Vulkan: pip install -r requirements-gpu.txt"
