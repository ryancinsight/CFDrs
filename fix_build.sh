#!/bin/bash

# Script to systematically fix build errors in cfd-core
source "/usr/local/cargo/env"

echo "Starting systematic repair of cfd-core..."

# Function to check if a file compiles
check_file() {
    local file=$1
    rustc --edition 2021 --crate-type lib "$file" 2>&1 | grep -q "error" && echo "ERROR: $file" || echo "OK: $file"
}

# Try to build and capture first error
echo "Attempting build to identify next error..."
cargo build --lib --package cfd-core 2>&1 | head -100

# List all files with errors
echo -e "\n=== Files with errors ==="
find /workspace/crates/cfd-core/src -name "*.rs" -exec sh -c 'rustc --edition 2021 --crate-type lib "$1" 2>&1 | grep -q "error" && basename "$1"' _ {} \; 2>/dev/null | sort

echo -e "\n=== Build attempt ==="
cargo build --workspace 2>&1 | grep "error:" | head -10