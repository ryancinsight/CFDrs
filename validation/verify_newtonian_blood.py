import math

D = 100e-6
L = 10e-3
dP = 4000
mu_newtonian = 0.0035  # Pa·s from source code

Q_analytical = (math.pi * D**4 * dP) / (128 * mu_newtonian * L)

print(f'Analytical Hagen-Poiseuille with μ=3.5 mPa·s:')
print(f'  Q = {Q_analytical*1e9:.6f} μL/s')
print(f'\nNumerical gave Q = 0.223883 μL/s')

error_pct = abs(Q_analytical - 0.223883e-9) / Q_analytical * 100
print(f'  Error: {error_pct:.4f}%')

if error_pct < 5:
    print('  ✓ PASS: Numerical matches analytical < 5%')
    print(f'\nCONCLUSION: "newtonian" blood = 3.5 mPa·s (NOT water!)')
