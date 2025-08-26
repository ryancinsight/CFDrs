#!/bin/bash

echo "Fixing import issues in all affected files..."

# Files that need fixing
FILES=(
    "crates/cfd-1d/src/resistance/models.rs"
    "crates/cfd-math/src/linear_solver/traits.rs"
    "crates/cfd-math/src/linear_solver/bicgstab.rs"
    "crates/cfd-math/src/linear_solver/conjugate_gradient.rs"
    "crates/cfd-math/src/linear_solver/preconditioners.rs"
)

for file in "${FILES[@]}"; do
    echo "Fixing $file..."
    
    # Add proper imports if not present
    if ! grep -q "use cfd_core::error::{Error, NumericalErrorKind}" "$file"; then
        # Find the first use statement and add our import after it
        sed -i '/^use /a use cfd_core::error::{Error, NumericalErrorKind};' "$file"
    fi
    
    # Fix the references
    sed -i 's/cfd_core::error::Error::Numerical(cfd_core::error::NumericalErrorKind/Error::Numerical(NumericalErrorKind/g' "$file"
done

echo "Import fixes applied."