#!/bin/bash
set -e

echo "=== Fixing ALL Copy Bounds Systematically ==="

# Fix all trait definitions
find crates -name "*.rs" -type f -exec sed -i \
    -e 's/trait \([A-Za-z]*\)<T: RealField>/trait \1<T: RealField + Copy>/g' \
    -e 's/trait \([A-Za-z]*\)<T: RealField + Send>/trait \1<T: RealField + Copy + Send>/g' \
    -e 's/trait \([A-Za-z]*\)<T: RealField + Sync>/trait \1<T: RealField + Copy + Sync>/g' \
    {} \;

# Fix all impl blocks
find crates -name "*.rs" -type f -exec sed -i \
    -e 's/impl<T: RealField>/impl<T: RealField + Copy>/g' \
    -e 's/impl<T: RealField + FromPrimitive>/impl<T: RealField + FromPrimitive + Copy>/g' \
    -e 's/impl<T: RealField + Debug>/impl<T: RealField + Debug + Copy>/g' \
    -e 's/impl<T: RealField + Display>/impl<T: RealField + Display + Copy>/g' \
    -e 's/impl<T: RealField + Float>/impl<T: RealField + Float + Copy>/g' \
    {} \;

# Fix struct definitions
find crates -name "*.rs" -type f -exec sed -i \
    -e 's/struct \([A-Za-z]*\)<T: RealField>/struct \1<T: RealField + Copy>/g' \
    -e 's/enum \([A-Za-z]*\)<T: RealField>/enum \1<T: RealField + Copy>/g' \
    {} \;

# Fix function signatures
find crates -name "*.rs" -type f -exec sed -i \
    -e 's/fn \([a-z_]*\)<T: RealField>/fn \1<T: RealField + Copy>/g' \
    -e 's/pub fn \([a-z_]*\)<T: RealField>/pub fn \1<T: RealField + Copy>/g' \
    {} \;

# Fix where clauses
find crates -name "*.rs" -type f -exec sed -i \
    -e 's/where T: RealField,/where T: RealField + Copy,/g' \
    -e 's/where T: RealField$/where T: RealField + Copy/g' \
    {} \;

# Prevent double Copy
find crates -name "*.rs" -type f -exec sed -i \
    -e 's/Copy + Copy/Copy/g' \
    -e 's/+ Copy + Copy/+ Copy/g' \
    {} \;

echo "=== Copy bounds fixed ==="