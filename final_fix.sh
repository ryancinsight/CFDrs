#!/bin/bash

echo "=== Final comprehensive fixes ==="

# Fix interpolation.rs line 299
sed -i '299s/acc \* (x - xj) \/ (self.x_data\[i\] - xj)/acc * (x - *xj) \/ (self.x_data[i] - *xj)/' crates/cfd-math/src/interpolation.rs

# Fix interpolation.rs line 310
sed -i '310s/yi \* self/\*yi * self/' crates/cfd-math/src/interpolation.rs

# Add clones where needed for matrix operations
sed -i 's/let error = (&y_next - &y_new)/let error = (\&y_next - \&y_new)/' crates/cfd-core/src/time.rs

# Build
cd /workspace
source /usr/local/cargo/env
cargo build --all 2>&1 | tail -20