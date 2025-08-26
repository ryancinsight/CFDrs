#!/bin/bash

# Comprehensive script to fix all structural issues in the codebase

source "/usr/local/cargo/env"

echo "Starting comprehensive repair of CFD Suite..."

# Function to fix a Rust file
fix_rust_file() {
    local file=$1
    echo "Fixing: $file"
    
    # Create a temporary backup
    cp "$file" "${file}.bak"
    
    # Read file content
    content=$(cat "$file")
    
    # Count braces
    open_count=$(echo "$content" | grep -o '{' | wc -l)
    close_count=$(echo "$content" | grep -o '}' | wc -l)
    
    # If more closing than opening braces, remove excess from end
    if [ $close_count -gt $open_count ]; then
        excess=$((close_count - open_count))
        # Remove trailing closing braces
        for i in $(seq 1 $excess); do
            content=$(echo "$content" | sed '$ {/^[[:space:]]*}[[:space:]]*$/d;}')
        done
        echo "$content" > "$file"
    fi
    
    # Try to compile and check for errors
    rustc --edition 2021 --crate-type lib "$file" 2>/dev/null
    if [ $? -eq 0 ]; then
        echo "  ✓ Fixed successfully"
        rm "${file}.bak"
    else
        echo "  ⚠ Needs manual intervention"
        mv "${file}.bak" "$file"
    fi
}

# Fix all Rust files in cfd-core
echo "Fixing cfd-core module..."
find /workspace/crates/cfd-core/src -name "*.rs" -type f | while read file; do
    fix_rust_file "$file"
done

# Fix all Rust files in other modules
for module in cfd-math cfd-mesh cfd-1d cfd-2d cfd-3d cfd-io cfd-validation; do
    echo "Fixing $module module..."
    find /workspace/crates/$module/src -name "*.rs" -type f 2>/dev/null | while read file; do
        fix_rust_file "$file"
    done
done

echo "Attempting full build..."
cargo build --workspace 2>&1 | tail -20

echo "Repair attempt complete."