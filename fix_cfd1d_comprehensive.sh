#!/bin/bash

echo "=== Comprehensive cfd-1d fixes ==="

# Add Copy bounds to all generic types
echo "Adding Copy bounds..."
find crates/cfd-1d -name "*.rs" -exec sed -i \
  -e 's/<T: RealField>/<T: RealField + Copy>/g' \
  -e 's/<T: RealField + From/<T: RealField + Copy + From/g' \
  -e 's/<T: RealField + Send/<T: RealField + Copy + Send/g' {} \;

# Fix arithmetic with references
echo "Fixing arithmetic operations..."
find crates/cfd-1d -name "*.rs" -exec sed -i \
  -e 's/\* viscosity/\* *viscosity/g' \
  -e 's/\* density/\* *density/g' \
  -e 's/\* diameter/\* *diameter/g' \
  -e 's/\* length/\* *length/g' \
  -e 's/\* width/\* *width/g' \
  -e 's/\* height/\* *height/g' \
  -e 's/\/ diameter/\/ *diameter/g' \
  -e 's/\/ area/\/ *area/g' {} \;

# Fix HashMap get operations
echo "Fixing HashMap operations..."
find crates/cfd-1d -name "*.rs" -exec sed -i \
  's/properties\.get(&id)/properties.get(\&id).cloned()/g' {} \;

# Fix network cloning
echo "Fixing network operations..."
sed -i 's/self\.network/self.network.clone()/g' crates/cfd-1d/src/channel.rs

echo "=== Done ==="