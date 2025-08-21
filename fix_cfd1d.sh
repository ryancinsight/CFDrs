#!/bin/bash
echo "Fixing cfd-1d issues..."

# Add Copy bounds
find crates/cfd-1d -name "*.rs" -exec sed -i 's/<T: RealField>/<T: RealField + Copy>/g' {} \;

# Fix arithmetic
find crates/cfd-1d -name "*.rs" -exec sed -i -e 's/\* viscosity/\* *viscosity/g' -e 's/\* density/\* *density/g' {} \;

echo "Done"
