#!/usr/bin/env python3
"""Test blood type fix in Poiseuille3DSolver"""

import pycfdrs
import math

D = 100e-6  # 100 μm
L = 10e-3   # 10 mm
dP = -4000  # Pa (negative for forward flow)

solver = pycfdrs.Poiseuille3DSolver(D, L, 20, 16, 30)

print("Testing blood type implementations:")
print("="*60)

# Test newtonian
try:
    r1 = solver.solve(dP, 'newtonian')
    print(f"✓ Newtonian: Q = {r1.flow_rate*1e12:.6f} pL/s")
except Exception as e:
    print(f"✗ Newtonian ERROR: {e}")

# Test casson
try:
    r2 = solver.solve(dP, 'casson')
    print(f"✓ Casson: Q = {r2.flow_rate*1e12:.6f} pL/s")
except Exception as e:
    print(f"✗ Casson ERROR: {e}")

# Test carreau_yasuda
try:
    r3 = solver.solve(dP, 'carreau_yasuda')
    print(f"✓ Carreau-Yasuda: Q = {r3.flow_rate*1e12:.6f} pL/s")
except Exception as e:
    print(f"✗ Carreau-Yasuda ERROR: {e}")

# Test invalid type
try:
    r4 = solver.solve(dP, 'invalid')
    print(f"✗ Should have errored for invalid type!")
except Exception as e:
    print(f"✓ Invalid type correctly rejected: {e}")

print("\n" + "="*60)
print("Expected values (based on μ at γ̇=100 s⁻¹):")
print("="*60)

mu_n = 0.0035    # Pa·s
mu_cas = 0.0043851  # Pa·s  
mu_car = 0.0047077  # Pa·s

Q_n = (math.pi * D**4 * abs(dP)) / (128 * mu_n * L)
Q_cas = (math.pi * D**4 * abs(dP)) / (128 * mu_cas * L)
Q_car = (math.pi * D**4 * abs(dP)) / (128 * mu_car * L)

print(f"Newtonian (μ=3.50 mPa·s): Q = {Q_n*1e12:.6f} pL/s")
print(f"Casson (μ=4.39 mPa·s): Q = {Q_cas*1e12:.6f} pL/s")
print(f"Carreau-Yasuda (μ=4.71 mPa·s): Q = {Q_car*1e12:.6f} pL/s")
