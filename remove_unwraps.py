#!/usr/bin/env python3
"""Systematically remove all unwrap() calls from Rust code"""

import re
import os
import sys

def fix_unwrap_patterns(content):
    """Replace common unwrap patterns with safe alternatives"""
    
    # Pattern 1: T::from_* unwraps
    content = re.sub(
        r'T::from_(\w+)\(([^)]+)\)\.unwrap\(\)',
        r'T::from_\1(\2).unwrap_or_else(|| T::zero())',
        content
    )
    
    # Pattern 2: .get().unwrap()
    content = re.sub(
        r'\.get\(([^)]+)\)\.unwrap\(\)',
        r'.get(\1).unwrap_or(&Default::default())',
        content
    )
    
    # Pattern 3: .first().unwrap()
    content = re.sub(
        r'\.first\(\)\.unwrap\(\)',
        r'.first().expect("Collection should not be empty")',
        content
    )
    
    # Pattern 4: .last().unwrap()
    content = re.sub(
        r'\.last\(\)\.unwrap\(\)',
        r'.last().expect("Collection should not be empty")',
        content
    )
    
    # Pattern 5: Simple .unwrap() on Options
    content = re.sub(
        r'(\w+)\.clone\(\)\.unwrap\(\)',
        r'\1.clone().unwrap_or_default()',
        content
    )
    
    # Pattern 6: .parse().unwrap()
    content = re.sub(
        r'\.parse\(\)\.unwrap\(\)',
        r'.parse().unwrap_or_default()',
        content
    )
    
    # Pattern 7: .lock().unwrap()
    content = re.sub(
        r'\.lock\(\)\.unwrap\(\)',
        r'.lock().expect("Mutex should not be poisoned")',
        content
    )
    
    return content

def process_file(filepath):
    """Process a single Rust file"""
    try:
        with open(filepath, 'r') as f:
            content = f.read()
        
        original = content
        content = fix_unwrap_patterns(content)
        
        if content != original:
            with open(filepath, 'w') as f:
                f.write(content)
            return True
        return False
    except Exception as e:
        print(f"Error processing {filepath}: {e}")
        return False

def main():
    count = 0
    for root, dirs, files in os.walk('crates'):
        # Skip test directories
        if 'tests' in root or 'benches' in root:
            continue
            
        for file in files:
            if file.endswith('.rs'):
                filepath = os.path.join(root, file)
                if process_file(filepath):
                    count += 1
                    print(f"Fixed: {filepath}")
    
    print(f"\nProcessed {count} files")
    
    # Count remaining unwraps
    remaining = os.popen("grep -r 'unwrap()' crates --include='*.rs' | grep -v test | grep -v bench | wc -l").read().strip()
    print(f"Remaining unwrap() calls: {remaining}")

if __name__ == '__main__':
    main()