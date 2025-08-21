#!/bin/bash

echo "=== Final comprehensive fixes ==="

# Fix cfd-2d remaining issues
echo "Fixing cfd-2d..."

# Fix fdm.rs
sed -i 's/rhs\[linear_idx\] = source_val/rhs[linear_idx] = source_val/' crates/cfd-2d/src/fdm.rs

# Fix fvm.rs arithmetic
sed -i 's/gradient \* dx/*gradient * dx/g' crates/cfd-2d/src/fvm.rs
sed -i 's/gradient \* dy/*gradient * dy/g' crates/cfd-2d/src/fvm.rs
sed -i 's/rhs\[linear_idx\] = source_value/rhs[linear_idx] = *source_value/' crates/cfd-2d/src/fvm.rs

# Fix lbm.rs
sed -i 's/= boundary_value;/= *boundary_value;/g' crates/cfd-2d/src/lbm.rs
sed -i 's/velocity \* dt/*velocity * dt/g' crates/cfd-2d/src/lbm.rs

# Fix cfd-1d remaining issues
echo "Fixing cfd-1d..."

# Add missing Copy bounds
find crates/cfd-1d -name "*.rs" -exec sed -i \
  's/<T: RealField>/<T: RealField + Copy>/g' {} \;

# Fix arithmetic with references
find crates/cfd-1d -name "*.rs" -exec sed -i \
  -e 's/\([a-z_][a-z0-9_]*\) \* \([a-z_][a-z0-9_]*\)/*\1 * *\2/g' \
  -e 's/\([a-z_][a-z0-9_]*\) \/ \([a-z_][a-z0-9_]*\)/*\1 \/ *\2/g' {} \;

# Fix HashMap operations
sed -i 's/properties\.get(&id)/properties.get(\&id).cloned()/g' crates/cfd-1d/src/network.rs

# Fix channel.rs
sed -i 's/self\.network/self.network.clone()/g' crates/cfd-1d/src/channel.rs

echo "=== Done ==="