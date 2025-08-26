#!/usr/bin/env python3
"""
Aggressive fix for corrupted Rust files.
This script identifies and fixes common corruption patterns.
"""

import os
import re
from pathlib import Path

def fix_rust_file_aggressive(filepath):
    """Aggressively fix structural issues in a Rust file."""
    try:
        with open(filepath, 'r') as f:
            lines = f.readlines()
        
        # Remove excessive closing braces at the end
        while lines and lines[-1].strip() == '}':
            lines.pop()
        
        content = ''.join(lines)
        
        # Count braces
        open_count = content.count('{')
        close_count = content.count('}')
        
        # Add missing closing braces
        if open_count > close_count:
            missing = open_count - close_count
            content += '\n' + '}\n' * missing
        
        # Fix common patterns
        # Pattern 1: impl blocks without closing brace
        content = re.sub(
            r'(impl[^{]*\{[^}]*?)(\nimpl\s+|\npub\s+struct\s+|\npub\s+enum\s+|\npub\s+trait\s+|\n#\[|\Z)',
            r'\1\n}\n\2',
            content
        )
        
        # Pattern 2: struct/enum definitions without closing brace
        content = re.sub(
            r'(pub\s+(?:struct|enum)[^{]*\{[^}]*?)(\npub\s+|\nimpl\s+|\n#\[|\Z)',
            r'\1\n}\n\2',
            content
        )
        
        # Pattern 3: functions without closing brace
        content = re.sub(
            r'(fn\s+\w+[^{]*\{[^}]*?)(\n\s*(?:pub\s+)?fn\s+|\nimpl\s+|\npub\s+(?:struct|enum)\s+|\n}\s*\n|\Z)',
            r'\1\n    }\n\2',
            content
        )
        
        # Write back
        with open(filepath, 'w') as f:
            f.write(content)
        
        return True
    except Exception as e:
        print(f"Error processing {filepath}: {e}")
        return False

def main():
    """Fix all Rust files in the workspace."""
    workspace = Path('/workspace/crates')
    fixed_count = 0
    
    for rust_file in workspace.rglob('*.rs'):
        if fix_rust_file_aggressive(rust_file):
            print(f"Fixed: {rust_file}")
            fixed_count += 1
    
    print(f"\nFixed {fixed_count} files")

if __name__ == "__main__":
    main()