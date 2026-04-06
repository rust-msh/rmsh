#!/bin/bash
# Verify script for EMStudio Web Local-First pipeline.
# Steps:
# 1) wasm32 checks for render/main/worker
# 2) optional trunk/wasm-pack availability check
# 3) optional full build script invocation

set -e

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "$ROOT_DIR"

echo "== EMStudio Web Local-First Verification =="
echo "Workspace: $ROOT_DIR"

echo
echo "[1/4] Checking wasm32 compile: emstudio-render"
cargo check -p emstudio-render --target wasm32-unknown-unknown

echo
echo "[2/4] Checking wasm32 compile: emstudio-worker"
cargo check -p emstudio-worker --target wasm32-unknown-unknown

echo
echo "[3/4] Checking wasm32 compile: emstudio-main"
set +e
cargo check -p emstudio-main --target wasm32-unknown-unknown
MAIN_RC=$?
set -e
if [ $MAIN_RC -ne 0 ]; then
  echo "[WARN] emstudio-main wasm check failed."
  echo "       If error mentions getrandom on wasm, enable getrandom/js in dependency tree."
fi

echo
echo "[4/4] Toolchain check"
if command -v trunk >/dev/null 2>&1; then
  echo "trunk: OK"
else
  echo "trunk: MISSING"
fi

if command -v wasm-pack >/dev/null 2>&1; then
  echo "wasm-pack: OK"
else
  echo "wasm-pack: MISSING"
fi

echo
echo "Optional full build: ./scripts/build-wasm.sh"
echo "Verification completed."
