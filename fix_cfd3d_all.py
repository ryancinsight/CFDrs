"""Comprehensive fix for all remaining from_f64, from(), and method ambiguity in cfd-3d."""
import re
import os

BASE = r"c:\Users\RyanClanton\gcli\software\millifluidic_design\CFDrs\crates\cfd-3d\src"

ALL_RS_FILES = []
for root, dirs, files in os.walk(BASE):
    for f in files:
        if f.endswith('.rs'):
            ALL_RS_FILES.append(os.path.join(root, f))

FP = '<T as FromPrimitive>'

def fix_file(path):
    with open(path, 'r', encoding='utf-8') as f:
        content = f.read()
    original = content

    # Fix T::from_f64( -> <T as FromPrimitive>::from_f64(
    if 'T::from_f64(' in content:
        content = content.replace('T::from_f64(', f'{FP}::from_f64(')
        # Make sure FromPrimitive is imported
        if 'FromPrimitive' not in content.split('use ')[0] if 'use ' in content else True:
            # It's likely already imported via num_traits::FromPrimitive or similar
            pass

    # Fix T::from(<f64_expr>) -> <T as From<f64>>::from(<f64_expr>)  
    content = re.sub(
        r'\bT::from\(([^)]+)\)',
        r'<T as From<f64>>::from(\1)',
        content
    )

    # Fix T::max_value() -> <T as nalgebra::RealField>::max_value()
    content = re.sub(r'\bT::max_value\(\)', '<T as nalgebra::RealField>::max_value()', content)
    # Fix T::min_value() -> <T as nalgebra::RealField>::min_value()
    content = re.sub(r'\bT::min_value\(\)', '<T as nalgebra::RealField>::min_value()', content)

    # Fix .abs() on T type -> Float::abs(x)
    # Pattern: (expr).abs() where we need to disambiguate
    # Only fix in generic contexts - replace with Float::abs(...)
    content = re.sub(
        r'([\w\d_\(\)]+)\.abs\(\)',
        r'num_traits::Float::abs(\1)',
        content
    )

    # Fix .sqrt() on T type
    content = re.sub(
        r'([\w\d_\(\)\[\]\.]+)\.sqrt\(\)',
        r'num_traits::Float::sqrt(\1)',
        content
    )

    # Fix .min(expr) on T type
    content = re.sub(
        r'(\w+)\.min\((\w+)\)',
        r'num_traits::Float::min(\1, \2)',
        content
    )

    # Fix .max(expr) on T type
    content = re.sub(
        r'(\w+)\.max\((\w+)\)',
        r'num_traits::Float::max(\1, \2)',
        content
    )

    # Add FromPrimitive import if we use it but it's not imported
    if f'{FP}' in content and 'use num_traits::FromPrimitive;' not in content and 'use num_traits::{Float, FromPrimitive' not in content and 'use num_traits::{FromPrimitive' not in content:
        # Check if we can find num_traits usage to inject after
        content = re.sub(
            r'(use num_traits::Float;)',
            r'use num_traits::{Float, FromPrimitive};',
            content,
            count=1
        )
        # If that pattern didn't match, try another
        if 'use num_traits::{Float, FromPrimitive}' not in content and f'{FP}' in content:
            content = re.sub(
                r'(use num_traits[^;]+;)',
                r'\1\nuse num_traits::FromPrimitive;',
                content,
                count=1
            )

    if content != original:
        with open(path, 'w', encoding='utf-8') as f:
            f.write(content)
        return True
    return False

changed = 0
for path in ALL_RS_FILES:
    if fix_file(path):
        changed += 1
        print(f"Fixed: {os.path.relpath(path, BASE)}")

print(f"\nTotal files changed: {changed}")
