import os

def fix_syntax(path):
    with open(path, 'r', encoding='utf-8') as f:
        content = f.read()
    orig = content
    
    content = content.replace('(T::one() + T::one()))', '(T::one() + T::one())')
    content = content.replace('(T::one() + T::one() + T::one()))', '(T::one() + T::one() + T::one())')
    content = content.replace('(T::one() + T::one() + T::one() + T::one()))', '(T::one() + T::one() + T::one() + T::one())')
    content = content.replace('T::pi())', 'T::pi()')
    content = content.replace('T::pi());', 'T::pi();')
    content = content.replace('expect("Mathematical constant conversion compromised"))', 'expect("Mathematical constant conversion compromised")')
    
    if content != orig:
        with open(path, 'w', encoding='utf-8') as f:
            f.write(content)
        print("Fixed syntax in", path)

rootDir = r"c:\Users\RyanClanton\gcli\software\millifluidic_design\CFDrs\crates\cfd-1d\src"
for subdir, dirs, files in os.walk(rootDir):
    for f in files:
        if f.endswith('.rs'):
            fix_syntax(os.path.join(subdir, f))
