#!/bin/bash
set -e

echo "=== Removing unnecessary clone() calls ==="

# Remove .clone() on numeric types (now that we have Copy)
find crates -name "*.rs" -type f -exec sed -i \
    -e 's/\.clone()//g' \
    {} \;

# Re-add clone only for String and complex types that actually need it
find crates -name "*.rs" -type f -exec sed -i \
    -e 's/\([a-z_]*\)\.id\([^a-z_]\)/\1.id.clone()\2/g' \
    -e 's/\([a-z_]*\)\.name\([^a-z_]\)/\1.name.clone()\2/g' \
    -e 's/String::from(\([^)]*\))/String::from(\1)/g' \
    {} \;

echo "=== Clone calls removed ==="