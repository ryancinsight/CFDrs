"""Fix specific method ambiguity in stabilization.rs, shape_functions.rs, ibm/interpolation.rs."""
import re
import os

BASE = r"c:\Users\RyanClanton\gcli\software\millifluidic_design\CFDrs\crates\cfd-3d\src"

def fix_method_ambiguities(content):
    """Fix ambiguous method calls for generic T types on specific patterns."""
    # Fix .min(expr) -> Float::min(x, expr) -- only when result is T typed
    content = re.sub(r'(\w+)\.min\((\w+)\)', r'num_traits::Float::min(\1, \2)', content)
    # Fix .max(expr) -> Float::max(x, expr)
    content = re.sub(r'(\w+)\.max\((\w+)\)', r'num_traits::Float::max(\1, \2)', content)
    # Fix (expr).max(T::zero()) 
    content = re.sub(r'\(([^)]+)\)\.max\(T::zero\(\)\)', r'num_traits::Float::max(\1, T::zero())', content)
    # Fix (expr).abs() -> Float::abs(expr)
    content = re.sub(r'\(([^)]+)\)\.abs\(\)', r'num_traits::Float::abs(\1)', content)
    # Fix sum.sqrt() -> Float::sqrt(sum)
    content = re.sub(r'(\w+)\.sqrt\(\)', r'num_traits::Float::sqrt(\1)', content)
    # Fix T::from(x) -> <T as From<f64>>::from(x)
    content = re.sub(r'\bT::from\(', r'<T as From<f64>>::from(', content)
    # Fix T::from_f64( -> <T as FromPrimitive>::from_f64(
    content = content.replace('T::from_f64(', '<T as FromPrimitive>::from_f64(')
    return content

def add_imports_if_needed(path, content):
    """Add missing imports if <T as FromPrimitive> is used but not imported."""
    needs_fp = '<T as FromPrimitive>' in content
    has_fp = 'FromPrimitive' in content.split('\n')[0:20]

    if needs_fp:
        lines = content.split('\n')
        has_fp_import = any('FromPrimitive' in l and l.strip().startswith('use') for l in lines[:30])
        if not has_fp_import:
            # Add after last use statement before first non-use line
            for i, line in enumerate(lines):
                if line.strip().startswith('use num_traits'):
                    if 'FromPrimitive' not in line:
                        if '{' in line and '}' in line:
                            # Replace "use num_traits::{A, B}" with "use num_traits::{A, B, FromPrimitive}"
                            content = content.replace(line, line.rstrip(';').rstrip('}') + ', FromPrimitive};', 1)
                        else:
                            # Add a new import line
                            insert_after = line
                            content = content.replace(insert_after, insert_after + '\nuse num_traits::FromPrimitive;', 1)
                    break
    return content

TARGETS = [
    r"fem\stabilization.rs",
    r"fem\shape_functions.rs",
    r"ibm\interpolation.rs",
    r"bifurcation\solver.rs",
    r"fem\quadrature.rs",
]

for rel in TARGETS:
    path = os.path.join(BASE, rel)
    if not os.path.exists(path):
        print(f"Skip: {rel}")
        continue
    with open(path, 'r', encoding='utf-8') as f:
        content = f.read()
    orig = content
    content = fix_method_ambiguities(content)
    content = add_imports_if_needed(path, content)
    if content != orig:
        with open(path, 'w', encoding='utf-8') as f:
            f.write(content)
        print(f"Fixed: {rel}")
    else:
        print(f"No changes: {rel}")
