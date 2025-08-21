#!/bin/bash

echo "=== Comprehensive cfd-2d fixes ==="

# Fix all reference dereferencing issues
echo "Fixing reference dereferencing..."
find crates/cfd-2d -name "*.rs" -exec sed -i \
  -e 's/= value;/= *value;/g' \
  -e 's/= gradient;/= *gradient;/g' \
  -e 's/= boundary_value;/= *boundary_value;/g' \
  -e 's/\([+-]\) gradient \*/\1 *gradient */g' \
  -e 's/\([+-]\) value \*/\1 *value */g' \
  -e 's/sum += value \*/sum += *value */g' \
  -e 's/diagonal = value/diagonal = *value/g' {} \;

# Fix DVector initialization
echo "Fixing DVector initialization..."
sed -i 's/DVector::zeros(n)/DVector::from_element(n, T::zero())/g' crates/cfd-2d/src/fdm.rs

# Fix HashMap get operations
echo "Fixing HashMap operations..."
sed -i 's/source\.get(&(i, j))/source.get(&(i, j)).copied()/g' crates/cfd-2d/src/fdm.rs
sed -i 's/boundary_values\.get(&(i, j))/boundary_values.get(&(i, j)).copied()/g' crates/cfd-2d/src/fdm.rs

# Fix arithmetic with references
echo "Fixing arithmetic operations..."
find crates/cfd-2d -name "*.rs" -exec sed -i \
  -e 's/\([a-z_][a-z0-9_]*\) \* dx/*\1 * dx/g' \
  -e 's/\([a-z_][a-z0-9_]*\) \* dy/*\1 * dy/g' \
  -e 's/dx \* \([a-z_][a-z0-9_]*\)/dx * *\1/g' \
  -e 's/dy \* \([a-z_][a-z0-9_]*\)/dy * *\1/g' {} \;

echo "=== Done ==="