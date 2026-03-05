import os, re

def clean_math(match):
    val_str = match.group(1)
    fallback = match.group(2)
    
    # Exact algebraic/trigonometric traits
    if val_str == 'PI' or val_str == 'std::f64::consts::PI': return 'T::pi()'
    if val_str == '1.0': return 'T::one()'
    if val_str == '0.0': return 'T::zero()'
    if val_str == '2.0': return '(T::one() + T::one())'
    if val_str == '3.0': return '(T::one() + T::one() + T::one())'
    if val_str == '4.0': return '(T::one() + T::one() + T::one() + T::one())'
    if val_str == '0.5': return '(T::one() / (T::one() + T::one()))'
    
    # Strict fallback unmasking
    return f'T::from_f64({val_str}).expect("Mathematical constant conversion compromised")'

rootDir = r"c:\Users\RyanClanton\gcli\software\millifluidic_design\CFDrs\crates\cfd-1d\src"
pattern = re.compile(r'T::from_f64\(([^)]+)\)\.(unwrap_or_else\([^)]+\)|unwrap_or\([^)]+\))')

for subdir, dirs, files in os.walk(rootDir):
    for f in files:
        if f.endswith('.rs'):
            path = os.path.join(subdir, f)
            with open(path, 'r', encoding='utf-8') as file: content = file.read()
            original = content
            content = pattern.sub(clean_math, content)
            if content != original:
                with open(path, 'w', encoding='utf-8') as file: file.write(content)
                print(f"Purified math in {path}")
