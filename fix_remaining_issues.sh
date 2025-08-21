#!/bin/bash

echo "=== Fixing remaining compilation issues ==="

# Add missing trait bounds to generic types
echo "Adding Copy bounds to generic types..."
find crates/cfd-1d crates/cfd-2d -name "*.rs" -exec sed -i \
  -e 's/<T: RealField>/<T: RealField + Copy>/g' \
  -e 's/<T: RealField + From/<T: RealField + Copy + From/g' \
  -e 's/<T: RealField + Send/<T: RealField + Copy + Send/g' {} \;

# Fix arithmetic operations with references
echo "Fixing arithmetic operations..."
find crates/cfd-1d crates/cfd-2d -name "*.rs" -exec sed -i \
  -e 's/\(&[a-z_][a-z0-9_]*\) \* T::/\*\1 * T::/g' \
  -e 's/T::.*() \* \(&[a-z_][a-z0-9_]*\)/T::from_f64(1.0).unwrap() * \*\2/g' {} \;

# Fix specific known issues
echo "Fixing specific issues..."

# Fix energy.rs
sed -i 's/temperature\[idx\] = boundary_temp;/temperature[idx] = *boundary_temp;/' crates/cfd-2d/src/energy.rs

# Fix fdm.rs
sed -i 's/T::zero()/T::zero()/' crates/cfd-2d/src/fdm.rs

# Fix channel.rs
sed -i 's/self\.network/self.network.clone()/' crates/cfd-1d/src/channel.rs

echo "=== Done ==="