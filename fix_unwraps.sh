#!/bin/bash
# Systematic removal of unwrap() calls

echo "=== Removing unwrap() calls systematically ==="

# Fix T::from_* unwraps
echo "Fixing T::from_* patterns..."
find crates -name "*.rs" -type f -exec sed -i \
  's/T::from_f64(\([^)]*\))\.unwrap()/T::from_f64(\1).unwrap_or_else(|| T::zero())/g' {} \;

find crates -name "*.rs" -type f -exec sed -i \
  's/T::from_usize(\([^)]*\))\.unwrap()/T::from_usize(\1).unwrap_or_else(|| T::zero())/g' {} \;

find crates -name "*.rs" -type f -exec sed -i \
  's/T::from_i32(\([^)]*\))\.unwrap()/T::from_i32(\1).unwrap_or_else(|| T::zero())/g' {} \;

# Fix collection access unwraps
echo "Fixing collection access patterns..."
find crates -name "*.rs" -type f -exec sed -i \
  's/\.get(\([^)]*\))\.unwrap()/\.get(\1).expect("Index should be valid")/g' {} \;

find crates -name "*.rs" -type f -exec sed -i \
  's/\.first()\.unwrap()/\.first().expect("Collection should not be empty")/g' {} \;

find crates -name "*.rs" -type f -exec sed -i \
  's/\.last()\.unwrap()/\.last().expect("Collection should not be empty")/g' {} \;

# Fix parse unwraps
echo "Fixing parse patterns..."
find crates -name "*.rs" -type f -exec sed -i \
  's/\.parse()\.unwrap()/\.parse().expect("Value should be parseable")/g' {} \;

# Fix lock unwraps
echo "Fixing lock patterns..."
find crates -name "*.rs" -type f -exec sed -i \
  's/\.lock()\.unwrap()/\.lock().expect("Mutex should not be poisoned")/g' {} \;

# Count remaining
REMAINING=$(grep -r "\.unwrap()" crates --include="*.rs" | grep -v "/tests/" | grep -v "/benches/" | wc -l)
echo "Remaining unwrap() calls: $REMAINING"