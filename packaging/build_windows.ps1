Param()
$ErrorActionPreference = 'Stop'

# Resolve repo root (script is in packaging/)
$ScriptDir = Split-Path -Parent $MyInvocation.MyCommand.Definition
$RepoRoot = Resolve-Path (Join-Path $ScriptDir "..")
Set-Location $RepoRoot

Write-Host "[packaging] Build: Windows"
python --version

python -m venv .venv-build
. .\.venv-build\Scripts\Activate.ps1

pip install --upgrade pip setuptools wheel
if (Test-Path -Path 'requirements.txt') {
  pip install -r requirements.txt
}
pip install pyinstaller

# cleanup
Remove-Item -Recurse -Force build, dist, *.spec -ErrorAction SilentlyContinue

pyinstaller --onefile --name genexplore main.py

New-Item -ItemType Directory -Path packaging\dist -Force | Out-Null
$exe = Join-Path $RepoRoot 'dist\genexplore.exe'
if (Test-Path $exe) {
  $zip = Join-Path $RepoRoot 'packaging\dist\genexplore-windows.zip'
  if (Test-Path $zip) { Remove-Item $zip -Force }
  Compress-Archive -Path $exe -DestinationPath $zip -Force
  Write-Host "[packaging] Artifact: $zip"
} else {
  Write-Host "[packaging] Build failed: dist\\genexplore.exe not found" -ForegroundColor Red
  Exit 1
}
