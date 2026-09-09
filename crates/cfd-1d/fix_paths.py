import os
import re

domain_modules = ['network', 'channel', 'components', 'junctions']
physics_modules = ['resistance', 'hemolysis', 'cell_separation', 'vascular']
solver_modules = ['analysis'] # core is handled specially

replacements = []
for mod in domain_modules:
    replacements.append((re.compile(rf'\bcrate::{mod}\b'), f'crate::domain::{mod}'))
    replacements.append((re.compile(rf'\bsuper::{mod}\b'), f'crate::domain::{mod}'))
    
for mod in physics_modules:
    replacements.append((re.compile(rf'\bcrate::{mod}\b'), f'crate::physics::{mod}'))
    replacements.append((re.compile(rf'\bsuper::{mod}\b'), f'crate::physics::{mod}'))

for mod in solver_modules:
    replacements.append((re.compile(rf'\bcrate::{mod}\b'), f'crate::solver::{mod}'))
    replacements.append((re.compile(rf'\bsuper::{mod}\b'), f'crate::solver::{mod}'))

# For solver, anything that was `crate::solver` and NOT `crate::solver::analysis` or `crate::solver::core`
replacements.append((re.compile(rf'\bcrate::solver(?!(::analysis|::core))\b'), 'crate::solver::core'))
replacements.append((re.compile(rf'\bsuper::solver(?!(::analysis|::core))\b'), 'crate::solver::core'))

rootDir = r"c:\Users\RyanClanton\gcli\software\millifluidic_design\CFDrs\crates\cfd-1d\src"

for subdir, dirs, files in os.walk(rootDir):
    for file in files:
        if file.endswith('.rs') and file != 'lib.rs':
            path = os.path.join(subdir, file)
            with open(path, 'r', encoding='utf-8') as f:
                content = f.read()

            new_content = content
            for (pattern, repl) in replacements:
                new_content = pattern.sub(repl, new_content)
                
            if new_content != content:
                with open(path, 'w', encoding='utf-8') as f:
                    f.write(new_content)
