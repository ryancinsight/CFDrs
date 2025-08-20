#!/bin/bash
# Script to replace unwrap() calls with proper error handling

echo "Replacing T::from_* unwrap() calls with ok_or_else..."

# Replace T::from_f64(x).unwrap() with proper error handling
find crates -name "*.rs" -type f -exec sed -i 's/T::from_f64(\([^)]*\))\.unwrap()/T::from_f64(\1).ok_or_else(|| cfd_core::error::Error::Numerical(cfd_core::error::NumericalErrorKind::InvalidFpOperation))?/g' {} \;

# Replace T::from_usize(x).unwrap() with proper error handling  
find crates -name "*.rs" -type f -exec sed -i 's/T::from_usize(\([^)]*\))\.unwrap()/T::from_usize(\1).ok_or_else(|| cfd_core::error::Error::Numerical(cfd_core::error::NumericalErrorKind::InvalidFpOperation))?/g' {} \;

# Replace T::from_i32(x).unwrap() with proper error handling
find crates -name "*.rs" -type f -exec sed -i 's/T::from_i32(\([^)]*\))\.unwrap()/T::from_i32(\1).ok_or_else(|| cfd_core::error::Error::Numerical(cfd_core::error::NumericalErrorKind::InvalidFpOperation))?/g' {} \;

echo "Replacing simple unwrap() calls on Option/Result..."

# Replace .get().unwrap() with .get().ok_or_else()
find crates -name "*.rs" -type f -exec sed -i 's/\.get(\([^)]*\))\.unwrap()/\.get(\1).ok_or_else(|| cfd_core::error::Error::InvalidConfiguration("Index out of bounds".into()))?/g' {} \;

# Replace .first().unwrap() with .first().ok_or_else()
find crates -name "*.rs" -type f -exec sed -i 's/\.first()\.unwrap()/\.first().ok_or_else(|| cfd_core::error::Error::InvalidConfiguration("Empty collection".into()))?/g' {} \;

# Replace .last().unwrap() with .last().ok_or_else()
find crates -name "*.rs" -type f -exec sed -i 's/\.last()\.unwrap()/\.last().ok_or_else(|| cfd_core::error::Error::InvalidConfiguration("Empty collection".into()))?/g' {} \;

echo "Done! Manual review still needed for complex unwrap() patterns."