#!/bin/bash
set -e

echo "=== Final comprehensive fix ==="

# Add Copy to EVERYTHING that uses RealField
find crates -name "*.rs" -type f -exec perl -i -pe \
    's/(<T:\s*RealField)([^>]*>)/$1 + Copy$2/g unless /\+ Copy/' \
    {} \;

# Fix specific patterns that were missed
find crates -name "*.rs" -type f -exec sed -i \
    -e 's/: RealField + Copy + Copy/: RealField + Copy/g' \
    -e 's/T: RealField,/T: RealField + Copy,/g' \
    -e 's/T: RealField>/T: RealField + Copy>/g' \
    -e 's/T: RealField + Copy + Copy/T: RealField + Copy/g' \
    {} \;

# Fix where clauses specifically
find crates -name "*.rs" -type f -exec sed -i \
    -e 's/where\s*T:\s*RealField\s*,/where T: RealField + Copy,/g' \
    -e 's/where\s*T:\s*RealField$/where T: RealField + Copy/g' \
    {} \;

# Fix impl blocks that might have been missed
find crates -name "*.rs" -type f -exec sed -i \
    -e 's/impl\s*<\s*T:\s*RealField\s*>/impl<T: RealField + Copy>/g' \
    {} \;

# Clean up any double Copy
find crates -name "*.rs" -type f -exec sed -i \
    -e 's/+ Copy + Copy/+ Copy/g' \
    -e 's/Copy + Copy/Copy/g' \
    {} \;

echo "=== Final fix complete ==="