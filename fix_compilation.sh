#!/bin/bash

# Fix compilation errors in the CFD codebase

echo "Fixing compilation errors..."

# Fix Copy bounds in numerical_methods.rs
sed -i 's/T: RealField + Clone/T: RealField + Copy/g' crates/cfd-core/src/domains/numerical_methods.rs

# Fix closure captures in numerical_methods.rs
sed -i 's/\.map(|val| val \* dt)/\.map(move |val| val * dt)/g' crates/cfd-core/src/domains/numerical_methods.rs

# Fix time.rs cloning issues
sed -i '181s/let mut y_new = state;/let mut y_new = state.clone();/' crates/cfd-core/src/time.rs
sed -i '182s/let y_old = state;/let y_old = state.clone();/' crates/cfd-core/src/time.rs
sed -i '195s/let y_next = &y_old + &(&f_val \* dt);/let y_next = y_old.clone() + f_val * dt;/' crates/cfd-core/src/time.rs

# Fix similar issues in Crank-Nicolson
sed -i '281s/let mut y_new = state;/let mut y_new = state.clone();/' crates/cfd-core/src/time.rs
sed -i '282s/let y_old = state;/let y_old = state.clone();/' crates/cfd-core/src/time.rs
sed -i '295s/let y_next = &y_old + &((&f_old + &f_new) \* half_dt);/let y_next = y_old.clone() + (f_old.clone() + f_new) * half_dt;/' crates/cfd-core/src/time.rs

# Fix factory.rs
sed -i '149s/self.factory.create(&typed_config)/self.factory.create(typed_config)/' crates/cfd-core/src/factory.rs

echo "Compilation fixes applied."