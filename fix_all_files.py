#!/usr/bin/env python3
"""
Script to systematically fix structural issues in Rust files.
Identifies and repairs common patterns from the corrupted refactoring.
"""

import os
import re
from pathlib import Path

def fix_rust_file(filepath):
    """Fix common structural issues in a Rust file."""
    try:
        with open(filepath, 'r') as f:
            content = f.read()
        
        original = content
        
        # Count braces
        open_braces = content.count('{')
        close_braces = content.count('}')
        
        # Remove excessive closing braces at end
        if close_braces > open_braces:
            lines = content.split('\n')
            # Remove lines that are just '}' at the end
            while lines and lines[-1].strip() == '}' and close_braces > open_braces:
                lines.pop()
                close_braces -= 1
            content = '\n'.join(lines)
        
        # Fix incomplete function bodies (pattern: missing closing brace)
        # Look for functions without proper closing
        content = re.sub(
            r'(\n\s*fn\s+\w+[^{]*\{[^}]*?)(\n\s*fn\s+|\n\s*pub\s+fn\s+|\n\s*impl\s+|\n\s*pub\s+struct\s+|\n\s*pub\s+enum\s+|\Z)',
            r'\1\n    }\n\2',
            content,
            flags=re.MULTILINE | re.DOTALL
        )
        
        # Fix incomplete impl blocks
        content = re.sub(
            r'(\n\s*impl[^{]*\{[^}]*?)(\n\s*impl\s+|\n\s*pub\s+struct\s+|\n\s*pub\s+enum\s+|\Z)',
            r'\1\n}\n\2',
            content,
            flags=re.MULTILINE | re.DOTALL
        )
        
        if content != original:
            with open(filepath, 'w') as f:
                f.write(content)
            return True
        return False
    except Exception as e:
        print(f"Error processing {filepath}: {e}")
        return False

def main():
    """Fix all Rust files in cfd-core."""
    fixed_count = 0
    cfd_core_path = Path('/workspace/crates/cfd-core/src')
    
    for rust_file in cfd_core_path.rglob('*.rs'):
        if fix_rust_file(rust_file):
            print(f"Fixed: {rust_file}")
            fixed_count += 1
    
    print(f"\nFixed {fixed_count} files")

if __name__ == "__main__":
    main()