#!/bin/bash
# WASM Build Script for EMStudio
# Builds both the main application and worker modules

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"
DIST_DIR="$PROJECT_ROOT/dist"

echo "🔨 EMStudio WASM Build"
echo "Project root: $PROJECT_ROOT"
echo "Output: $DIST_DIR"

# Clean distribution directory
if [ -d "$DIST_DIR" ]; then
    echo "🧹 Cleaning $DIST_DIR..."
    rm -rf "$DIST_DIR"
fi

mkdir -p "$DIST_DIR"

# Phase 1: Build worker WASM module
echo ""
echo "📦 Building worker WASM module..."
cd "$PROJECT_ROOT/crates/worker"

wasm-pack build \
    --target web \
    --out-dir "$DIST_DIR/worker" \
    --release

echo "✅ Worker WASM built to $DIST_DIR/worker"

# Phase 2: Build main application via trunk
echo ""
echo "🎯 Building main application with trunk..."
cd "$PROJECT_ROOT/crates/main"

trunk build --release --public-url /

echo "✅ Main app built to $PROJECT_ROOT/crates/main/dist"

# Phase 3: Copy worker files to main dist (if needed)
# Trunk's post-build hooks will handle this, but we can do it manually too
if [ -d "$PROJECT_ROOT/crates/main/dist" ]; then
    MAIN_DIST="$PROJECT_ROOT/crates/main/dist"
    echo ""
    echo "📂 Organizing output..."
    
    # Create worker subdirectory if it doesn't exist
    mkdir -p "$MAIN_DIST/worker"
    
    # Copy worker WASM module
    if [ -f "$DIST_DIR/worker/emstudio_worker.js" ]; then
        cp "$DIST_DIR/worker/emstudio_worker.js" "$MAIN_DIST/worker/"
        cp "$DIST_DIR/worker/emstudio_worker_bg.wasm" "$MAIN_DIST/worker/"
        echo "✅ Worker module copied to $MAIN_DIST/worker"
    fi
    
    echo ""
    echo "📊 Build output:"
    echo "  HTML:                    $MAIN_DIST/index.html"
    echo "  Main WASM:              $MAIN_DIST/emstudio_main_bg.wasm"
    echo "  Main JS:                $MAIN_DIST/emstudio_main.js"
    echo "  Worker WASM:            $MAIN_DIST/worker/emstudio_worker_bg.wasm"
    echo "  Worker JS:              $MAIN_DIST/worker/emstudio_worker.js"
fi

echo ""
echo "✨ WASM build complete!"
echo ""
echo "📝 Next steps:"
echo "  1. Run a local web server: python3 -m http.server 8000 --directory $PROJECT_ROOT/crates/main/dist"
echo "  2. Open http://localhost:8000 in your browser"
echo "  3. Check browser console (F12) for any errors"
