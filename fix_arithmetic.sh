#!/bin/bash
set -e

echo "=== Fixing arithmetic operations with references ==="

# Fix all arithmetic operations in numerical_methods.rs
cd crates/cfd-core/src/domains

# Fix patterns where references need dereferencing
sed -i 's/\([a-z_]*\) \* dt/*\1 * dt/g' numerical_methods.rs
sed -i 's/\([a-z_]*\) + \([a-z_]*\) \* /*\1 + *\2 */g' numerical_methods.rs
sed -i 's/\([a-z_]*\) - \([a-z_]*\) \* /*\1 - *\2 */g' numerical_methods.rs
sed -i 's/\([a-z_]*\) \/ \([a-z_]*\)/*\1 \/ *\2/g' numerical_methods.rs

# Fix iterator patterns
sed -i 's/\.map(|x| x \* /\.map(|x| *x * /g' numerical_methods.rs
sed -i 's/\.map(|x| x + /\.map(|x| *x + /g' numerical_methods.rs
sed -i 's/\.map(|x| x - /\.map(|x| *x - /g' numerical_methods.rs
sed -i 's/\.map(|x| x \/ /\.map(|x| *x \/ /g' numerical_methods.rs

# Fix tuple patterns
sed -i 's/|(a, b)| a /|(a, b)| *a /g' numerical_methods.rs
sed -i 's/|(x, y)| x /|(x, y)| *x /g' numerical_methods.rs
sed -i 's/|(u, v)| u /|(u, v)| *u /g' numerical_methods.rs

# Clean up double dereferences
sed -i 's/\*\*/***/g' numerical_methods.rs
sed -i 's/\*T::/T::/g' numerical_methods.rs

cd /workspace

echo "=== Arithmetic fixes complete ==="