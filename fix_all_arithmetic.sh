#!/bin/bash

echo "=== Fixing all arithmetic operations systematically ==="

# Fix interpolation.rs line 299
echo "Fixing interpolation.rs..."
sed -i '299s/(x - \*xj)/(x - *xj)/' crates/cfd-math/src/interpolation.rs
sed -i '299s/(self\.x_data\[i\] - \*xj)/(self.x_data[i] - *xj)/' crates/cfd-math/src/interpolation.rs

# Fix linear_solver.rs line 210
echo "Fixing linear_solver.rs..."
sed -i '210s/sum = sum + val \* z\[\*j\]/sum = sum + *val * z[*j]/' crates/cfd-math/src/linear_solver.rs

# Fix sparse.rs line 215
echo "Fixing sparse.rs..."
sed -i '215s/value \* x\[col_idx\]/*value * x[col_idx]/' crates/cfd-math/src/sparse.rs

# Fix sparse.rs line 224
echo "Fixing sparse.rs line 224..."
sed -i '224s/\.map(|v| v \* v)/\.map(|v| *v * *v)/' crates/cfd-math/src/sparse.rs

# Fix all similar patterns
echo "Applying general fixes..."
find crates/cfd-math -name "*.rs" -exec sed -i \
  -e 's/\([a-z_][a-z0-9_]*\)\.\([a-z_][a-z0-9_]*\) \* \([a-z_][a-z0-9_]*\)\[\([^]]*\)\]/\1.\2 * \3[\4]/' \
  -e 's/sum = sum + \([a-z_][a-z0-9_]*\) \*/sum = sum + *\1 */' {} \; 2>/dev/null || true

echo "=== Fixes applied ==="