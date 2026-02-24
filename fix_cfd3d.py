import os

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
        
        with open(filepath, 'w', encoding='utf-8') as f:
            f.write(content)
    except Exception as e:
        print(f"Error processing {filepath}: {e}")

for root, _, files in os.walk('crates/cfd-3d/src'):
    for file in files:
        if file.endswith('.rs'):
            fix_file(os.path.join(root, file))

print("Fixed imports in cfd-3d")
