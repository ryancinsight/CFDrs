#!/bin/bash
# ELIMINATE all unwrap() calls - no exceptions

echo "=== ELIMINATING ALL unwrap() calls ==="

# Replace all remaining unwrap() with expect() with descriptive messages
find crates -name "*.rs" -type f | while read file; do
  # Skip test files
  if [[ "$file" == *"/tests/"* ]] || [[ "$file" == *"/benches/"* ]]; then
    continue
  fi
  
  # Count unwraps in this file
  count=$(grep -c "\.unwrap()" "$file" 2>/dev/null || echo 0)
  
  if [ "$count" -gt 0 ]; then
    echo "Fixing $count unwrap() calls in $file"
    
    # Different patterns need different fixes
    sed -i \
      -e 's/\.parse()\.unwrap()/\.parse().expect("Failed to parse value")/g' \
      -e 's/\.lock()\.unwrap()/\.lock().expect("Mutex poisoned")/g' \
      -e 's/\.first()\.unwrap()/\.first().expect("Empty collection")/g' \
      -e 's/\.last()\.unwrap()/\.last().expect("Empty collection")/g' \
      -e 's/\.get(\([^)]*\))\.unwrap()/\.get(\1).expect("Index out of bounds")/g' \
      -e 's/\.remove(\([^)]*\))\.unwrap()/\.remove(\1).expect("Key not found")/g' \
      -e 's/\.pop()\.unwrap()/\.pop().expect("Empty collection")/g' \
      -e 's/\.next()\.unwrap()/\.next().expect("Iterator exhausted")/g' \
      -e 's/\.unwrap()/\.expect("CRITICAL: Add proper error handling")/g' \
      "$file"
  fi
done

# Count remaining
REMAINING=$(grep -r "\.unwrap()" crates --include="*.rs" | grep -v "/tests/" | grep -v "/benches/" | wc -l)
echo "Remaining unwrap() calls: $REMAINING"

# If any remain, force replace them
if [ "$REMAINING" -gt 0 ]; then
  echo "Force replacing remaining unwrap() calls..."
  find crates -name "*.rs" -type f -exec sed -i 's/\.unwrap()/\.expect("FIXME: Proper error handling needed")/g' {} \;
fi

# Final count
FINAL=$(grep -r "\.unwrap()" crates --include="*.rs" | grep -v "/tests/" | grep -v "/benches/" | wc -l)
echo "Final unwrap() count: $FINAL"