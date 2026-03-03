import os
import re

def process_directory(directory):
    count = 0
    pattern_f64 = re.compile(r'(T::from_(?:f64|f32|u32|i32)(?:\([^)]+\)))\.unwrap\(\)')
    
    for root, _, files in os.walk(directory):
        for file in files:
            if file.endswith('.rs'):
                filepath = os.path.join(root, file)
                with open(filepath, 'r', encoding='utf-8') as f:
                    content = f.read()
                
                # Replace T::from_f64(x).unwrap() with T::from_f64(x).unwrap_or_else(num_traits::Zero::zero)
                new_content = pattern_f64.sub(r'\1.unwrap_or_else(num_traits::Zero::zero)', content)
                
                if new_content != content:
                    with open(filepath, 'w', encoding='utf-8') as f:
                        f.write(new_content)
                    print(f"Updated {filepath}")
                    count += 1
    return count

if __name__ == "__main__":
    crates_dir = os.path.join("crates")
    c = process_directory(crates_dir)
    print(f"Total files updated: {c}")
