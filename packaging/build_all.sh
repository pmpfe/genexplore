#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

case "$(uname -s)" in
  Linux*) bash "$SCRIPT_DIR/build_linux.sh" ;;
  Darwin*) bash "$SCRIPT_DIR/build_macos.sh" ;;
  MINGW*|MSYS*|CYGWIN*) pwsh -File "$SCRIPT_DIR/build_windows.ps1" ;;
  *) echo "Unsupported OS: $(uname -s)"; exit 1 ;;
esac
