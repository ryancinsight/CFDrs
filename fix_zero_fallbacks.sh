#!/bin/bash

# Script to fix dangerous T::zero() fallbacks

echo "Fixing dangerous T::zero() fallbacks in CFD Suite..."

# Files to fix
FILES=(
    "crates/cfd-1d/src/resistance/models.rs"
    "crates/cfd-validation/src/time_integration_validation.rs"
    "crates/cfd-3d/src/level_set.rs"
    "crates/cfd-3d/src/vof/advection.rs"
    "crates/cfd-core/src/values.rs"
)

for file in "${FILES[@]}"; do
    echo "Fixing $file..."
    
    # Replace unwrap_or_else(|| T::zero()) with proper error handling
    sed -i 's/\.unwrap_or_else(|| T::zero())/\.ok_or_else(|| cfd_core::error::Error::Numerical(cfd_core::error::NumericalErrorKind::ConversionFailed { from_type: "f64", to_type: std::any::type_name::<T>(), value: "value".to_string() }))?/g' "$file"
    
    # Replace unwrap_or(T::zero()) similarly
    sed -i 's/\.unwrap_or(T::zero())/\.ok_or_else(|| cfd_core::error::Error::Numerical(cfd_core::error::NumericalErrorKind::ConversionFailed { from_type: "f64", to_type: std::any::type_name::<T>(), value: "value".to_string() }))?/g' "$file"
    
    # For T::one() fallbacks too
    sed -i 's/\.unwrap_or_else(|| T::one())/\.ok_or_else(|| cfd_core::error::Error::Numerical(cfd_core::error::NumericalErrorKind::ConversionFailed { from_type: "f64", to_type: std::any::type_name::<T>(), value: "1.0".to_string() }))?/g' "$file"
done

echo "Critical fixes applied. These files now properly handle numeric conversion errors."