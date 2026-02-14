#!/usr/bin/env python3
"""
Fix Unicode encoding issues in validation scripts for Windows compatibility
"""

import os
import glob

# Unicode character replacements
replacements = {
    '\u03bc': 'um',  # μ (micro)
    '\u00d7': 'x',    # × (multiplication)
    '\u2713': '[OK]', # ✓ (checkmark)
    '\u274c': '[FAIL]', # ❌ (cross mark)
}

def fix_file(filepath):
    """Fix Unicode characters in a single file"""
    with open(filepath, 'r', encoding='utf-8') as f:
        content = f.read()
    
    # Check if file needs fixing
    needs_fix = any(char in content for char in replacements.keys())
    
    if not needs_fix:
        return False
    
    # Apply replacements
    for old, new in replacements.items():
        content = content.replace(old, new)
    
    # Add UTF-8 encoding handling at the top if not already present
    if 'sys.stdout = io.TextIOWrapper' not in content:
        lines = content.split('\n')
        insert_pos = 0
        
        # Find position after imports
        for i, line in enumerate(lines):
            if line.startswith('import ') or line.startswith('from '):
                insert_pos = i + 1
        
        # Insert encoding handling
        encoding_code = [
            '',
            '# Set UTF-8 encoding for stdout to handle Unicode characters',
            'import io',
            'if sys.platform == "win32":',
            '    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8")',
            '    sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding="utf-8")',
        ]
        
        lines = lines[:insert_pos] + encoding_code + lines[insert_pos:]
        content = '\n'.join(lines)
    
    # Write fixed content
    with open(filepath, 'w', encoding='utf-8') as f:
        f.write(content)
    
    return True

def main():
    """Fix all validation Python scripts"""
    script_dir = os.path.dirname(__file__)
    pattern = os.path.join(script_dir, 'complete_*.py')
    
    fixed_count = 0
    for filepath in glob.glob(pattern):
        filename = os.path.basename(filepath)
        if fix_file(filepath):
            print(f"Fixed: {filename}")
            fixed_count += 1
        else:
            print(f"Skipped (no Unicode): {filename}")
    
    print(f"\nFixed {fixed_count} files")

if __name__ == '__main__':
    main()
