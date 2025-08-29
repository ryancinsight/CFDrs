#!/bin/bash
# Migration script to replace deprecated Fluid with ConstantPropertyFluid

echo "Migrating deprecated Fluid to ConstantPropertyFluid..."

# Replace Fluid imports
find crates -name "*.rs" -type f -exec sed -i \
  -e 's/use crate::fluid::Fluid;/use crate::fluid::ConstantPropertyFluid;/g' \
  -e 's/use cfd_core::fluid::Fluid;/use cfd_core::fluid::ConstantPropertyFluid;/g' \
  {} \;

# Replace Fluid type annotations  
find crates -name "*.rs" -type f -exec sed -i \
  -e 's/: Fluid</: ConstantPropertyFluid</g' \
  -e 's/-> Fluid</-> ConstantPropertyFluid</g' \
  -e 's/Fluid::/ConstantPropertyFluid::/g' \
  {} \;

# Replace .create() with .new()
find crates -name "*.rs" -type f -exec sed -i \
  's/ConstantPropertyFluid::create/ConstantPropertyFluid::new/g' {} \;

echo "Migration complete. Please review changes and run cargo build."