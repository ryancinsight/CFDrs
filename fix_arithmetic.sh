#!/bin/bash

echo "=== Fixing arithmetic operations ==="

# Fix interpolation.rs
echo "Fixing interpolation.rs..."
sed -i '112s/(x - x0) \/ (x1 - x0)/(x - *x0) \/ (*x1 - *x0)/' crates/cfd-math/src/interpolation.rs
sed -i '113s/y0 + t \* (y1 - y0)/*y0 + t * (*y1 - *y0)/' crates/cfd-math/src/interpolation.rs

# Fix all reference arithmetic issues
find crates -name "*.rs" -exec sed -i \
  -e 's/\(&[a-zA-Z_][a-zA-Z0-9_]*\) - \(&[a-zA-Z_][a-zA-Z0-9_]*\)/(*\1 - *\2)/g' \
  -e 's/\(&[a-zA-Z_][a-zA-Z0-9_]*\) + \([a-zA-Z_][a-zA-Z0-9_]*\)/*\1 + \2/g' \
  -e 's/\(&[a-zA-Z_][a-zA-Z0-9_]*\) \* \([a-zA-Z_][a-zA-Z0-9_]*\)/*\1 * \2/g' {} \; 2>/dev/null || true

echo "=== Arithmetic fixes applied ==="