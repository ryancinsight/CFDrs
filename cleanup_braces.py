#!/usr/bin/env python3
"""
Clean up extra closing braces that were added by the aggressive fix.
"""

import re
from pathlib import Path

def cleanup_extra_braces(filepath):
    """Remove extra closing braces after struct/enum definitions."""
    try:
        with open(filepath, 'r') as f:
            content = f.read()
        
        # Pattern: struct/enum with fields ending with }}\n}
        # This indicates an extra closing brace
        content = re.sub(
            r'(pub\s+(?:struct|enum)[^{]*\{[^}]*?\n\s*\}\n)\s*\}\n',
            r'\1',
            content
        )
        
        # Pattern: Remove double }} at end of impl blocks
        content = re.sub(
            r'(\n\s*\}\n)\s*\}\n\}',
            r'\1}',
            content
        )
        
        with open(filepath, 'w') as f:
            f.write(content)
        
        return True
    except Exception as e:
        print(f"Error processing {filepath}: {e}")
        return False

def main():
    """Clean up all Rust files in cfd-core."""
    workspace = Path('/workspace/crates/cfd-core/src')
    fixed_count = 0
    
    for rust_file in workspace.rglob('*.rs'):
        if cleanup_extra_braces(rust_file):
            fixed_count += 1
    
    print(f"Cleaned {fixed_count} files")

if __name__ == "__main__":
    main()