#!/bin/bash

echo "=== Fixing all compilation errors ==="

# Fix missing Copy bounds in cfd-math
echo "Fixing cfd-math Copy bounds..."
find crates/cfd-math -name "*.rs" -exec sed -i \
  -e 's/T: RealField + FromPrimitive/T: RealField + FromPrimitive + Copy/g' \
  -e 's/T: RealField + From/T: RealField + From<f64> + Copy/g' \
  -e 's/T: RealField,/T: RealField + Copy,/g' \
  -e 's/T: RealField>/T: RealField + Copy>/g' {} \;

# Fix missing Copy bounds in cfd-core
echo "Fixing cfd-core Copy bounds..."
find crates/cfd-core -name "*.rs" -exec sed -i \
  -e 's/T: RealField + FromPrimitive/T: RealField + FromPrimitive + Copy/g' \
  -e 's/T: RealField,/T: RealField + Copy,/g' \
  -e 's/T: RealField>/T: RealField + Copy>/g' {} \;

# Fix factory.rs Config type
echo "Fixing factory.rs..."
sed -i '19s/type Config: SolverConfiguration<T>;/type Config: SolverConfiguration<T> + Clone;/' crates/cfd-core/src/factory.rs
sed -i '149s/create(typed_config)/create(typed_config.clone())/' crates/cfd-core/src/factory.rs

# Fix all other crates
echo "Fixing remaining crates..."
for crate in cfd-1d cfd-2d cfd-3d cfd-mesh cfd-io cfd-validation; do
  echo "  Fixing $crate..."
  find crates/$crate -name "*.rs" -exec sed -i \
    -e 's/T: RealField + FromPrimitive/T: RealField + FromPrimitive + Copy/g' \
    -e 's/T: RealField,/T: RealField + Copy,/g' \
    -e 's/T: RealField>/T: RealField + Copy>/g' {} \;
done

echo "=== Compilation fixes applied ==="