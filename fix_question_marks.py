#!/usr/bin/env python3
"""Fix ? operators in functions that don't return Result"""

import re
import os
import sys

def fix_file(filepath):
    with open(filepath, 'r') as f:
        lines = f.readlines()
    
    modified = False
    in_function = False
    returns_result = False
    function_indent = 0
    
    for i in range(len(lines)):
        line = lines[i]
        
        # Detect function signatures
        fn_match = re.match(r'^(\s*)(pub\s+)?(fn\s+\w+.*?)\s*->\s*(.*?)(\s*\{)?$', line)
        if fn_match:
            function_indent = len(fn_match.group(1))
            return_type = fn_match.group(4)
            returns_result = 'Result' in return_type or 'Option' in return_type
            in_function = True
            continue
        
        # Check if we're still in the function
        if in_function:
            indent = len(line) - len(line.lstrip())
            if indent <= function_indent and line.strip():
                in_function = False
                returns_result = False
                continue
        
        # Fix ? operators in non-Result functions
        if in_function and not returns_result and '?' in line:
            # Replace .ok_or_else(...)? with .unwrap_or_else(...)
            new_line = re.sub(
                r'\.ok_or_else\([^)]+\)\?',
                lambda m: '.unwrap_or_else(' + m.group(0)[13:-2] + ')',
                line
            )
            
            # For T::from_* patterns, use unwrap_or with default
            new_line = re.sub(
                r'T::from_(\w+)\(([^)]+)\)\.ok_or_else\([^)]+\)\?',
                r'T::from_\1(\2).unwrap_or(T::zero())',
                new_line
            )
            
            if new_line != line:
                lines[i] = new_line
                modified = True
    
    if modified:
        with open(filepath, 'w') as f:
            f.writelines(lines)
        return True
    return False

def main():
    count = 0
    for root, dirs, files in os.walk('crates'):
        for file in files:
            if file.endswith('.rs'):
                filepath = os.path.join(root, file)
                if fix_file(filepath):
                    count += 1
                    print(f"Fixed: {filepath}")
    
    print(f"\nFixed {count} files")

if __name__ == '__main__':
    main()