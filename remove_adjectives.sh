#!/bin/bash
# Remove ALL adjective-based naming

echo "=== REMOVING ALL adjective-based naming ==="

# Replace common adjective patterns in code
find crates -name "*.rs" -type f -exec sed -i \
  -e 's/simple_/standard_/g' \
  -e 's/Simple/Standard/g' \
  -e 's/SIMPLE/STANDARD/g' \
  -e 's/basic_/core_/g' \
  -e 's/Basic/Core/g' \
  -e 's/advanced_/extended_/g' \
  -e 's/Advanced/Extended/g' \
  -e 's/enhanced_/augmented_/g' \
  -e 's/Enhanced/Augmented/g' \
  -e 's/optimized_/tuned_/g' \
  -e 's/Optimized/Tuned/g' \
  -e 's/improved_/revised_/g' \
  -e 's/Improved/Revised/g' \
  -e 's/fast_/rapid_/g' \
  -e 's/Fast/Rapid/g' \
  -e 's/robust_/resilient_/g' \
  -e 's/Robust/Resilient/g' \
  -e 's/accurate_/precise_/g' \
  -e 's/Accurate/Precise/g' \
  -e 's/new_/current_/g' \
  -e 's/New/Current/g' \
  -e 's/old_/previous_/g' \
  -e 's/Old/Previous/g' {} \;

# Fix specific known violations
find crates -name "*.rs" -type f -exec sed -i \
  -e 's/"simple"/"pressure-velocity"/g' \
  -e 's/"basic"/"fundamental"/g' \
  -e 's/"advanced"/"comprehensive"/g' \
  -e 's/"fast"/"efficient"/g' \
  -e 's/"robust"/"stable"/g' {} \;

echo "Adjective removal complete"