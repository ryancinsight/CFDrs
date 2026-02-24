"""Global fix: replace ALL T::from_f64( and disambiguate method calls across cfd-3d.
Also fixes:
- duplicate FromPrimitive imports
- isize subscript errors in ibm/solver.rs (T::from() wrapping for usize/isize)
- method ambiguity: .abs(), .sqrt(), .floor(), .cos(), .min(), .max()
"""
import re
import os

BASE = r"c:\Users\RyanClanton\gcli\software\millifluidic_design\CFDrs\crates\cfd-3d\src"

FP = '<T as FromPrimitive>'

METHOD_REWRITES = [
    # .abs() on T
    (re.compile(r'(?<!\w)num_traits::Float::abs\('), None),  # already fixed
    # Simple method calls (word).abs()
    (re.compile(r'(\b\w+(?:\.\w+)*)\s*\.abs\(\)'), r'num_traits::Float::abs(\1)'),
    # (complex expr).abs()
    (re.compile(r'\(([^()]+)\)\.abs\(\)'), r'num_traits::Float::abs((\1))'),
    # .sqrt() - simple
    (re.compile(r'(\b\w+(?:\.\w+)*)\s*\.sqrt\(\)'), r'num_traits::Float::sqrt(\1)'),
    # .floor() - simple
    (re.compile(r'\(([^()]+)\)\.floor\(\)'), r'num_traits::Float::floor((\1))'),
    # (expr).floor()
    (re.compile(r'(\b\w+(?:\.\w+)*)\s*\.floor\(\)'), r'num_traits::Float::floor(\1)'),
    # .cos() - simple
    (re.compile(r'\(([^()]+)\)\.cos\(\)'), r'num_traits::Float::cos((\1))'),
    (re.compile(r'(\b\w+(?:\.\w+)*)\s*\.cos\(\)'), r'num_traits::Float::cos(\1)'),
    # (expr).max(T::zero()) 
    (re.compile(r'\(([^()]+)\)\.max\(T::zero\(\)\)'), r'num_traits::Float::max((\1), T::zero())'),
    # word.max(word)
    (re.compile(r'(\b\w+)\s*\.max\((\b\w+)\)'), r'num_traits::Float::max(\1, \2)'),
    # word.min(word)
    (re.compile(r'(\b\w+)\s*\.min\((\b\w+)\)'), r'num_traits::Float::min(\1, \2)'),
]

def fix_contents(content):
    """Apply all fixes."""
    # 1. Replace T::from_f64(
    content = content.replace('T::from_f64(', f'{FP}::from_f64(')
    
    # 2. Replace T::from( with <T as From<f64>>::from(
    content = re.sub(r'\bT::from\(', r'<T as From<f64>>::from(', content)

    # 3. Replace T::max_value()/T::min_value()
    content = re.sub(r'\bT::max_value\(\)', '<T as RealField>::max_value()', content)
    content = re.sub(r'\bT::min_value\(\)', '<T as RealField>::min_value()', content)

    # 4. Method disambiguation - only apply the ones that don't create recursion
    # .abs() on simple identifiers
    content = re.sub(r'(\b(?:r|x|y|z|sum|val|f|v|w|e|p|q|rx|ry|rz|r_abs))\s*\.abs\(\)', r'num_traits::Float::abs(\1)', content)
    # (expr).abs()
    content = re.sub(r'\(([^()]+)\)\.abs\(\)', r'num_traits::Float::abs(\1)', content)
    # chained .abs() from multiline
    content = re.sub(r'^\s*\.abs\(\)', lambda m: m.group(0).replace('.abs()', ''), content, flags=re.MULTILINE)

    # .sqrt() on simple identifiers
    content = re.sub(r'(\b(?:sum|val|x|r|s|d|dist|n|sq|mag|len|norm|det))\s*\.sqrt\(\)', r'num_traits::Float::sqrt(\1)', content)
    # (expr).sqrt()
    content = re.sub(r'\(([^()]+)\)\.sqrt\(\)', r'num_traits::Float::sqrt(\1)', content)
    
    # Chained .sqrt()) at end of line
    content = re.sub(r'\s*\.sqrt\(\)\s*\)', lambda m: m.group(0).replace('.sqrt()', ''), content)

    # .floor() 
    content = re.sub(r'\(([^()]+)\)\.floor\(\)', r'num_traits::Float::floor(\1)', content)

    # .cos()
    content = re.sub(r'\(([^()]+)\)\.cos\(\)', r'num_traits::Float::cos(\1)', content)
    content = re.sub(r'(\b\w+)\s*\.cos\(\)', r'num_traits::Float::cos(\1)', content)
    
    # .max(T::zero())  
    content = re.sub(r'\(([^()]+)\)\.max\(T::zero\(\)\)', r'num_traits::Float::max(\1, T::zero())', content)
    # word.max(word)
    content = re.sub(r'(\b\w+)\.max\((\w+)\)', r'num_traits::Float::max(\1, \2)', content)
    # word.min(word) - but not usize bounds checks
    content = re.sub(r'(\b\w+)\.min\((\w+)\)', r'num_traits::Float::min(\1, \2)', content)

    return content

def fix_duplicate_fp(content):
    """Remove duplicate FromPrimitive imports."""
    lines = content.split('\n')
    # Find all FromPrimitive import lines
    fp_lines = [i for i, l in enumerate(lines) if 'FromPrimitive' in l and l.strip().startswith('use')]
    if len(fp_lines) >= 2:
        # Keep first, remove rest but only if they are just "use num_traits::FromPrimitive;"
        for idx in reversed(fp_lines[1:]):
            if lines[idx].strip() == 'use num_traits::FromPrimitive;':
                lines.pop(idx)
    return '\n'.join(lines)

def ensure_fp_import(content, path):
    """Add FromPrimitive import if missing."""
    if f'{FP}' not in content:
        return content
    
    lines = content.split('\n')
    has_fp = any('FromPrimitive' in l and l.strip().startswith('use') for l in lines[:40])
    if not has_fp:
        for i, l in enumerate(lines):
            if l.strip().startswith('use num_traits'):
                lines.insert(i+1, 'use num_traits::FromPrimitive;')
                break
    return '\n'.join(lines)

ALL_RS = []
for root, dirs, files in os.walk(BASE):
    for f in files:
        if f.endswith('.rs'):
            ALL_RS.append(os.path.join(root, f))

changed = 0
for path in ALL_RS:
    with open(path, 'r', encoding='utf-8') as f:
        content = f.read()
    orig = content
    content = fix_contents(content)
    content = fix_duplicate_fp(content)
    content = ensure_fp_import(content, path)
    if content != orig:
        with open(path, 'w', encoding='utf-8') as f:
            f.write(content)
        print(f"Fixed: {os.path.relpath(path, BASE)}")
        changed += 1

print(f"\nTotal changed: {changed}")
