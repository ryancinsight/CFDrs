#!/usr/bin/env python3
"""Check what Poiseuille classes are available in pycfdrs"""

import pycfdrs

print("Available Poiseuille-related attributes in pycfdrs:")
poiseuille_attrs = [attr for attr in dir(pycfdrs) if 'poiseuille' in attr.lower()]
for attr in sorted(poiseuille_attrs):
    obj = getattr(pycfdrs, attr)
    print(f"  {attr}: {type(obj)}")

print("\nAll attributes:")
all_attrs = [attr for attr in dir(pycfdrs) if not attr.startswith('_')]
for attr in sorted(all_attrs):
    print(f"  {attr}")
