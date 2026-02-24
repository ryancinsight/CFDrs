"""Fix corrupted struct field accesses (resolution.num_traits::Float::max(...)) 
   and do proper disambiguations via direct string search/replace."""
import re
import os

BASE = r"c:\Users\RyanClanton\gcli\software\millifluidic_design\CFDrs\crates\cfd-3d\src"

def revert_corruption(content):
    """Revert cases where regex incorrectly replaced field.method() with field.num_traits::Float::method()."""
    # This pattern: word.num_traits::Float::method( should be reverted to just word.method(
    # BUT only when it's a struct field access, not a standalone Float call
    # Detect: something like self.config.resolution.num_traits::Float::max(0, 1)
    # vs the correct usage: num_traits::Float::max(a, b) 
    
    # Fix: word.num_traits::Float::max(x, y) -> word.max(x, y)  (was corrupted)
    content = re.sub(
        r'(\w+)\.num_traits::Float::max\(([^)]+)\)',
        r'\1.max(\2)',
        content
    )
    content = re.sub(
        r'(\w+)\.num_traits::Float::min\(([^)]+)\)',
        r'\1.min(\2)',
        content
    )
    content = re.sub(
        r'(\w+)\.num_traits::Float::abs\(',
        r'\1.abs(',
        content
    )
    content = re.sub(
        r'(\w+)\.num_traits::Float::sqrt\(',
        r'\1.sqrt(',
        content
    )
    content = re.sub(
        r'(\w+)\.num_traits::Float::floor\(',
        r'\1.floor(',
        content
    )
    content = re.sub(
        r'(\w+)\.num_traits::Float::cos\(',
        r'\1.cos(',
        content
    )
    content = re.sub(
        r'(\w+)\.num_traits::Float::sin\(',
        r'\1.sin(',
        content
    )
    content = re.sub(
        r'(\w+)\.num_traits::Float::powi\(',
        r'\1.powi(',
        content
    )
    content = re.sub(
        r'(\w+)\.num_traits::Float::powf\(',
        r'\1.powf(',
        content
    )
    return content

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
    content = revert_corruption(content)
    if content != orig:
        with open(path, 'w', encoding='utf-8') as f:
            f.write(content)
        print(f"Reverted: {os.path.relpath(path, BASE)}")
        changed += 1

print(f"Total reverted: {changed}")
