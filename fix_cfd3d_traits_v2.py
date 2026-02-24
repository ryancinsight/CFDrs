import os
import re

def fix_file(filepath):
    try:
        with open(filepath, 'r', encoding='utf-8') as f:
            content = f.read()

        # Reverse previous bad trait injections
        content = content.replace(' + cfd_mesh::domain::core::scalar::Scalar>', '>')
        content = content.replace('+ cfd_mesh::domain::core::scalar::Scalar>', '>')
        content = content.replace(' + cfd_mesh::domain::core::scalar::Scalar>(', '>(')
        
        # Inject cfd_mesh::domain::core::Scalar + at the start of ALL <T: ...>
        # Match <T: and replace with <T: cfd_mesh::domain::core::Scalar + 
        content = re.sub(r'<T:\s+', r'<T: cfd_mesh::domain::core::Scalar + ', content)
        content = content.replace('Float::epsilon()', 'num_traits::Float::epsilon()')

        with open(filepath, 'w', encoding='utf-8') as f:
            f.write(content)
    except Exception as e:
        print(f"Error processing {filepath}: {e}")

for root, _, files in os.walk('crates/cfd-3d/src'):
    for file in files:
        if file.endswith('.rs'):
            fix_file(os.path.join(root, file))

print("Applied trait bounds to cfd-3d v2")
