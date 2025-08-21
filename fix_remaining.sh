#!/bin/bash

echo "=== Fixing remaining compilation errors ==="

# Fix sparse.rs line 282
sed -i '282s/\.map(|(idx, val)| \*idx as f64 \* val)/\.map(|(idx, val)| *idx as f64 * *val)/' crates/cfd-math/src/sparse.rs

# Fix sparse.rs line 292
sed -i '292s/(val - mean)/(val - mean)/' crates/cfd-math/src/sparse.rs
sed -i '292s/\*val - mean/*val - mean/' crates/cfd-math/src/sparse.rs

# Fix differentiation.rs line 243
sed -i '243s/T::zero()/T::zero()/' crates/cfd-math/src/differentiation.rs

# Fix integration.rs line 157-158
sed -i '157s/f(a)/f(a)/' crates/cfd-math/src/integration.rs
sed -i '158s/val \* h/*val * h/' crates/cfd-math/src/integration.rs

# Fix vectorization.rs
sed -i '27s/a + b/*a + *b/' crates/cfd-math/src/vectorization.rs
sed -i '46s/a \* b/*a * *b/' crates/cfd-math/src/vectorization.rs
sed -i '65s/\.max()/\.max()/' crates/cfd-math/src/vectorization.rs
sed -i '90s/acc + val/acc + *val/' crates/cfd-math/src/vectorization.rs

echo "=== Fixes applied ==="