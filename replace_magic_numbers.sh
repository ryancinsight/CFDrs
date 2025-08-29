#!/bin/bash

# Replace magic number tolerances with named constants
# This script systematically replaces hardcoded numerical values with
# properly named constants from cfd_core::constants::solver_numerics

echo "Replacing magic numbers with named constants..."

# Replace 1e-6 (default tolerance)
find crates -name "*.rs" -type f | while read file; do
    # Skip test files and constants definitions
    if [[ ! "$file" =~ "test" ]] && [[ ! "$file" =~ "constants" ]]; then
        sed -i 's/T::from_f64(1e-6)/T::from_f64(cfd_core::constants::solver_numerics::DEFAULT_SOLVER_TOLERANCE)/g' "$file"
        sed -i 's/1e-6/cfd_core::constants::solver_numerics::DEFAULT_SOLVER_TOLERANCE/g' "$file"
    fi
done

# Replace 1e-8 (pressure tolerance)
find crates -name "*.rs" -type f | while read file; do
    if [[ ! "$file" =~ "test" ]] && [[ ! "$file" =~ "constants" ]]; then
        sed -i 's/T::from_f64(1e-8)/T::from_f64(cfd_core::constants::solver_numerics::PRESSURE_SOLVER_TOLERANCE)/g' "$file"
        sed -i 's/1e-8/cfd_core::constants::solver_numerics::PRESSURE_SOLVER_TOLERANCE/g' "$file"
    fi
done

# Replace 1e-10 (epsilon tolerance)
find crates -name "*.rs" -type f | while read file; do
    if [[ ! "$file" =~ "test" ]] && [[ ! "$file" =~ "constants" ]]; then
        sed -i 's/T::from_f64(1e-10)/T::from_f64(cfd_core::constants::solver_numerics::EPSILON_TOLERANCE)/g' "$file"
        sed -i 's/1e-10/cfd_core::constants::solver_numerics::EPSILON_TOLERANCE/g' "$file"
    fi
done

# Replace 1000 (max iterations)
find crates -name "*.rs" -type f | while read file; do
    if [[ ! "$file" =~ "test" ]] && [[ ! "$file" =~ "constants" ]]; then
        # Only replace standalone 1000, not as part of other numbers
        sed -i 's/max_iterations: 1000/max_iterations: cfd_core::constants::solver_numerics::DEFAULT_MAX_ITERATIONS/g' "$file"
    fi
done

echo "Magic number replacement complete"