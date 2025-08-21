#!/bin/bash

echo "=== Comprehensive Fix Script ==="

# Fix all reference arithmetic in interpolation.rs
echo "Fixing interpolation.rs..."
sed -i 's/\([a-z_][a-z0-9_]*\) - \([a-z_][a-z0-9_]*\)/(*\1 - *\2)/g' crates/cfd-math/src/interpolation.rs
sed -i 's/\([a-z_][a-z0-9_]*\) + \([a-z_][a-z0-9_]*\)/(*\1 + *\2)/g' crates/cfd-math/src/interpolation.rs
sed -i 's/\([a-z_][a-z0-9_]*\) \* \([a-z_][a-z0-9_]*\)/(*\1 * *\2)/g' crates/cfd-math/src/interpolation.rs
sed -i 's/\([a-z_][a-z0-9_]*\) \/ \([a-z_][a-z0-9_]*\)/(*\1 \/ *\2)/g' crates/cfd-math/src/interpolation.rs

# Fix dereferencing issues
sed -i 's/\*\*\([a-z_][a-z0-9_]*\)/*\1/g' crates/cfd-math/src/interpolation.rs

# Build and check
source /usr/local/cargo/env
cargo build --all 2>&1 | tail -20