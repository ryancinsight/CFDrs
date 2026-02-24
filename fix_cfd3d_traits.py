import os
import re

def fix_file(filepath):
    try:
        with open(filepath, 'r', encoding='utf-8') as f:
            content = f.read()

        # Fix imports
        content = content.replace('cfd_mesh::geometry::', 'cfd_mesh::domain::geometry::')
        content = content.replace('cfd_mesh::hierarchy::', 'cfd_mesh::application::hierarchy::')
        content = content.replace('cfd_mesh::topology::', 'cfd_mesh::domain::topology::')
        content = content.replace('use cfd_mesh::mesh::Mesh;', 'use cfd_mesh::IndexedMesh;')
        content = content.replace('cfd_mesh::mesh::Mesh', 'cfd_mesh::IndexedMesh')
        content = content.replace('mesh::Mesh<T>', 'IndexedMesh<T>')
        content = content.replace('&Mesh<T>', '&IndexedMesh<T>')
        
        # Any stray Mesh<T> used safely assuming we already took care of other meshes
        content = content.replace(' Mesh<T>', ' IndexedMesh<T>')
        content = content.replace('(mesh: IndexedMesh<T>', '(mesh: cfd_mesh::IndexedMesh<T>')
        content = content.replace('<Mesh<T>>', '<IndexedMesh<T>>')
        
        # Inject cfd_mesh::domain::core::scalar::Scalar trait bounds
        # Match `Float + From<f64>>` or similar
        content = re.sub(r'Float \+ From<f64>>', r'Float + From<f64> + cfd_mesh::domain::core::scalar::Scalar>', content)
        content = re.sub(r'Float> \(', r'Float + cfd_mesh::domain::core::scalar::Scalar> (', content)
        content = re.sub(r'Float>\(cell', r'Float + cfd_mesh::domain::core::scalar::Scalar>(cell', content)
        content = re.sub(r'Float>\(corners', r'Float + cfd_mesh::domain::core::scalar::Scalar>(corners', content)
        content = re.sub(r'Float>\(mesh', r'Float + cfd_mesh::domain::core::scalar::Scalar>(mesh', content)
        content = re.sub(r'Float::epsilon\(\)', r'num_traits::Float::epsilon()', content)

        with open(filepath, 'w', encoding='utf-8') as f:
            f.write(content)
    except Exception as e:
        print(f"Error processing {filepath}: {e}")

for root, _, files in os.walk('crates/cfd-3d/src'):
    for file in files:
        if file.endswith('.rs'):
            fix_file(os.path.join(root, file))

print("Applied trait bounds to cfd-3d")
