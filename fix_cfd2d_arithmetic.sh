#!/bin/bash

echo "Fixing cfd-2d arithmetic operations..."

# Fix energy.rs arithmetic operations
sed -i 's/\([a-z_][a-z0-9_]*\) \* \([a-z_][a-z0-9_]*\)\.\([a-z_][a-z0-9_]*\)/\*\1 * \2.\3/g' crates/cfd-2d/src/energy.rs
sed -i 's/\([a-z_][a-z0-9_]*\)\.\([a-z_][a-z0-9_]*\) \* \([a-z_][a-z0-9_]*\)/\1.\2 * \*\3/g' crates/cfd-2d/src/energy.rs

# Fix fdm.rs arithmetic operations
sed -i 's/val \* dx/*val * dx/g' crates/cfd-2d/src/fdm.rs
sed -i 's/val \* dy/*val * dy/g' crates/cfd-2d/src/fdm.rs

# Fix reference/value mismatches
find crates/cfd-2d -name "*.rs" -exec sed -i \
  -e 's/\(&[a-z_][a-z0-9_]*\) \* \([A-Z]\)/\*\1 * \2/g' \
  -e 's/\([A-Z]\) \* \(&[a-z_][a-z0-9_]*\)/\1 * \*\2/g' {} \;

echo "Done"