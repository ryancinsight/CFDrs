"""Fix all remaining from_f64 ambiguity and vertices.iter() tuple access across cfd-3d."""
import re
import os

BASE = r"c:\Users\RyanClanton\gcli\software\millifluidic_design\CFDrs\crates\cfd-3d\src"

# Files with from_f64 issues
FROM_F64_FILES = [
    r"fem\stabilization.rs",
    r"fem\quadrature.rs",
    r"fem\stress.rs",
    r"fem\solver.rs",
    r"fem\projection_solver.rs",
    r"bifurcation\solver.rs",
]

# Files with vertices.iter() -> (VertexId, &VertexData) tuples
VERTICES_ITER_FILES = [
    r"fem\solver.rs",
    r"fem\projection_solver.rs",
]

def fix_from_f64(content):
    # Replace T::from_f64( with <T as FromPrimitive>::from_f64(  
    content = content.replace('T::from_f64(', '<T as FromPrimitive>::from_f64(')
    # Replace T::max_value() with <T as RealField>::max_value()
    content = re.sub(r'T::max_value\(\)', '<T as RealField>::max_value()', content)
    # Replace T::min_value() with <T as RealField>::min_value()
    content = re.sub(r'T::min_value\(\)', '<T as RealField>::min_value()', content)
    # Fix .max(T::zero()) ambiguity -> use RealField::max
    content = re.sub(
        r'\(([^)]+)\)\.max\(T::zero\(\)\)',
        r'RealField::max(\1, T::zero())',
        content
    )
    return content

def fix_vertices_iter(content):
    # vertices.iter() returns (VertexId, &VertexData), not just &VertexData
    # Fix |v| v.position.coords -> |v| v.1.position.coords
    content = re.sub(
        r'\|v\| v\.position\.coords',
        '|v| v.1.position.coords',
        content
    )
    # Fix: let p = v.position.coords -> let p = v.1.position.coords
    content = re.sub(
        r'let p = v\.position\.coords;',
        'let p = v.1.position.coords;',
        content
    )
    return content

for rel in FROM_F64_FILES:
    path = os.path.join(BASE, rel)
    if not os.path.exists(path):
        print(f"Skip (not found): {rel}")
        continue
    with open(path, "r", encoding="utf-8") as f:
        content = f.read()
    new_content = fix_from_f64(content)
    if new_content != content:
        with open(path, "w", encoding="utf-8") as f:
            f.write(new_content)
        print(f"Fixed from_f64: {rel}")
    else:
        print(f"No from_f64 changes: {rel}")

for rel in VERTICES_ITER_FILES:
    path = os.path.join(BASE, rel)
    if not os.path.exists(path):
        print(f"Skip (not found): {rel}")
        continue
    with open(path, "r", encoding="utf-8") as f:
        content = f.read()
    new_content = fix_vertices_iter(content)
    if new_content != content:
        with open(path, "w", encoding="utf-8") as f:
            f.write(new_content)
        print(f"Fixed vertices.iter: {rel}")
    else:
        print(f"No vertices.iter changes: {rel}")

print("Done.")
