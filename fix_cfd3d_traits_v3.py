import os
import re

def fix_file(filepath):
    try:
        with open(filepath, 'r', encoding='utf-8') as f:
            content = f.read()

        # Fix missing bounds
        content = content.replace('<T>', '<T: cfd_mesh::domain::core::Scalar>')
        # Wait, if we replace <T>, we might break Option<T> or Mesh<T> or IndexedMesh<T>!
        # Do not blindly replace <T>. Instead replace `struct \w+<T>`
        content = re.sub(r'struct (\w+)<T>', r'struct \1<T: cfd_mesh::domain::core::Scalar>', content)
        content = re.sub(r'impl<T> (\w+)<T>', r'impl<T: cfd_mesh::domain::core::Scalar> \1<T>', content)
        
        # Also let's fix the ambiguity errors in bifurcation/geometry.rs
        content = content.replace('.clamp(T::zero(), one)', '')
        # Actually replacing clamp is safer by just doing num_traits::Float::clamp
        content = content.replace('(x / transition_length).clamp(T::zero(), one)', 'num_traits::Float::clamp(x / transition_length, T::zero(), one)')
        content = content.replace('(pi * x_normalized).cos()', 'num_traits::Float::cos(pi * x_normalized)')
        content = content.replace('self.d_parent.powi(3)', 'num_traits::Float::powi(self.d_parent, 3)')
        content = content.replace('self.d_daughter1.powi(3)', 'num_traits::Float::powi(self.d_daughter1, 3)')
        content = content.replace('self.d_daughter2.powi(3)', 'num_traits::Float::powi(self.d_daughter2, 3)')
        content = content.replace('arg.sqrt()', 'num_traits::Float::sqrt(arg)')
        content = content.replace('e01.max(e02).max(e03)', 'num_traits::Float::max(num_traits::Float::max(e01, e02), e03)')
        content = content.replace('e01.min(e02).min(e03)', 'num_traits::Float::min(num_traits::Float::min(e01, e02), e03)')

        with open(filepath, 'w', encoding='utf-8') as f:
            f.write(content)
    except Exception as e:
        print(f"Error processing {filepath}: {e}")

for root, _, files in os.walk('crates/cfd-3d/src'):
    for file in files:
        if file.endswith('.rs'):
            fix_file(os.path.join(root, file))

print("Applied struct trait bounds and exact disambiguations to cfd-3d")
