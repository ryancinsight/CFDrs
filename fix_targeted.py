"""Targeted fix: only specific disambiguation patterns reported by cargo check."""
import re
import os

BASE = r"c:\Users\RyanClanton\gcli\software\millifluidic_design\CFDrs\crates\cfd-3d\src"

FP = '<T as FromPrimitive>'

def fix_from_f64_and_from(path):
    with open(path, 'r', encoding='utf-8') as f:
        content = f.read()
    original = content

    # Fix T::from_f64( -> <T as FromPrimitive>::from_f64(
    content = content.replace('T::from_f64(', f'{FP}::from_f64(')

    # Fix T::from(<f64_expr>) -> <T as From<f64>>::from(<f64_expr>)
    content = re.sub(r'\bT::from\(', r'<T as From<f64>>::from(', content)

    # Ensure FromPrimitive is imported if <T as FromPrimitive> present
    if f'{FP}' in content:
        if 'use num_traits::FromPrimitive;' not in content \
           and '{Float, FromPrimitive}' not in content \
           and '{FromPrimitive,' not in content:
            # Try to extend existing num_traits Float import
            content = re.sub(
                r'use num_traits::Float;',
                'use num_traits::{Float, FromPrimitive};',
                content,
                count=1
            )
            if 'use num_traits::FromPrimitive;' not in content and '{FromPrimitive' not in content:
                # Try to add after first use num_traits line
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

# Files known to have these issues from the cargo output
PROBLEM_FILES = [
    r"fem\stress.rs",
    r"fem\stabilization.rs",
    r"fem\shape_functions.rs",
    r"fem\quadrature.rs",
    r"ibm\interpolation.rs",
    r"bifurcation\solver.rs",
    r"bifurcation\validation.rs",
    r"venturi\validation.rs",
    r"trifurcation\validation.rs",
    r"serpentine\validation.rs",
]

for rel in PROBLEM_FILES:
    path = os.path.join(BASE, rel)
    if not os.path.exists(path):
        print(f"Skip: {rel}")
        continue
    if fix_from_f64_and_from(path):
        print(f"Fixed: {rel}")
    else:
        print(f"No changes: {rel}")

print("Done.")
