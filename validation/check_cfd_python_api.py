#!/usr/bin/env python3
"""Check what Poiseuille classes are available in cfd_python"""

import cfd_python

print("Available Poiseuille-related attributes in cfd_python:")
poiseuille_attrs = [attr for attr in dir(cfd_python) if 'poiseuille' in attr.lower()]
for attr in sorted(poiseuille_attrs):
    obj = getattr(cfd_python, attr)
    print(f"  {attr}: {type(obj)}")

print("\nAll attributes:")
all_attrs = [attr for attr in dir(cfd_python) if not attr.startswith('_')]
for attr in sorted(all_attrs):
    print(f"  {attr}")
