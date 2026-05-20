#!/usr/bin/env bash
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$REPO_ROOT"

echo "[packaging] Build: Linux"
PYTHON=${PYTHON:-python3}
$PYTHON --version || true

# create venv
$PYTHON -m venv .venv-build
source .venv-build/bin/activate

pip install --upgrade pip setuptools wheel
if [ -f requirements.txt ]; then
  pip install -r requirements.txt
fi
pip install pyinstaller

# cleanup previous builds
rm -rf build dist __pycache__ *.spec

pyinstaller --onefile --name genexplore main.py

mkdir -p packaging/dist
tar -czf packaging/dist/genexplore-linux.tar.gz -C dist genexplore

echo "[packaging] Artifact: packaging/dist/genexplore-linux.tar.gz"
