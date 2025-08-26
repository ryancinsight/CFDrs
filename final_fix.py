#!/usr/bin/env python3
"""
Final targeted fix for remaining structural issues.
"""

import re
from pathlib import Path

def count_braces(content):
    """Count opening and closing braces."""
    opens = content.count('{')
    closes = content.count('}')
    return opens, closes

def fix_file_braces(filepath):
    """Fix brace mismatches in a file."""
    try:
        with open(filepath, 'r') as f:
            content = f.read()
        
        opens, closes = count_braces(content)
        
        if opens == closes:
            return False  # Already balanced
        
        if closes > opens:
            # Remove extra closing braces from the end
            lines = content.split('\n')
            while closes > opens and lines:
                if lines[-1].strip() == '}':
                    lines.pop()
                    closes -= 1
                else:
                    break
            content = '\n'.join(lines)
        elif opens > closes:
            # Add missing closing braces at the end
            missing = opens - closes
            content += '\n' + '}\n' * missing
        
        with open(filepath, 'w') as f:
            f.write(content)
        
        return True
    except Exception as e:
        print(f"Error processing {filepath}: {e}")
        return False

def main():
    """Fix all Rust files with brace issues."""
    files_to_fix = [
        '/workspace/crates/cfd-core/src/aggregates.rs',
        '/workspace/crates/cfd-core/src/domains/fluid_dynamics/rans.rs',
    ]
    
    # Also check all cfd-core files
    workspace = Path('/workspace/crates/cfd-core/src')
    for rust_file in workspace.rglob('*.rs'):
        files_to_fix.append(str(rust_file))
    
    fixed_count = 0
    for filepath in set(files_to_fix):  # Remove duplicates
        if fix_file_braces(filepath):
            print(f"Fixed: {filepath}")
            fixed_count += 1
    
    print(f"\nFixed {fixed_count} files")

if __name__ == "__main__":
    main()