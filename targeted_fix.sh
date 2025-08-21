#!/bin/bash
set -e
echo "=== Applying targeted fixes ==="
# Fix all RealField without Copy
find crates -name "*.rs" -type f -exec sed -i \
    -e 's/\(<T: RealField\)\([^+>]*>\)/\1 + Copy\2/g' \
    -e 's/+ Copy + Copy/+ Copy/g' \
    {} \;
# Fix where clauses
find crates -name "*.rs" -type f -exec sed -i \
    -e 's/where T: RealField\([^+]\)/where T: RealField + Copy\1/g' \
    -e 's/where T: RealField$/where T: RealField + Copy/g' \
    {} \;
echo "=== Targeted fixes applied ==="
