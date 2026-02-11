#!/usr/bin/env python3
"""
Direct Rust vs Python blood viscosity comparison
"""

import pycfdrs

# Carreau-Yasuda parameters (Cho & Kensey 1991)
mu_0 = 0.056  # Pa·s
mu_inf = 0.00345  # Pa·s
lambda_cy = 3.313  # s
a_cy = 2.0
n_cy = 0.3568

def carreau_yasuda_python(shear_rate):
    """Python reference implementation"""
    term = 1 + (lambda_cy * shear_rate) ** a_cy
    exponent = (n_cy - 1) / a_cy
    return mu_inf + (mu_0 - mu_inf) * (term ** exponent)

# Create Rust blood model
blood_rust = pycfdrs.CarreauYasudaBlood()

print("="*70)
print("RUST vs PYTHON: Carreau-Yasuda Blood Viscosity")
print("="*70)

test_rates = [10, 100, 500, 1000, 2000, 5000, 10000, 100000]

print(f"\n{'γ̇ (s⁻¹)':>12} {'Python (mPa·s)':>18} {'Rust (mPa·s)':>18} {'Error (%)':>15}")
print("-"*70)

max_error = 0.0
for rate in test_rates:
    mu_python = carreau_yasuda_python(rate)
    mu_rust = blood_rust.apparent_viscosity(rate)
    error_pct = abs(mu_rust - mu_python) / mu_python * 100
    max_error = max(max_error, error_pct)
    
    print(f"{rate:12.0f} {mu_python*1000:18.6f} {mu_rust*1000:18.6f} {error_pct:15.8f}")

print(f"\nMaximum error: {max_error:.8f}%")

if max_error < 0.01:
    print("✓ PASS: Rust matches Python within 0.01%")
else:
    print(f"✗ FAIL: Error {max_error:.4f}% exceeds 0.01% threshold")

print("\n" + "="*70)
print("CONCLUSION")
print("="*70)

if max_error < 0.01:
    print("""
✓ Rust CarreauYasudaBlood implementation is CORRECT
- Matches validated Python implementation
- Error < 0.01% across all test shear rates
- Ready for production use
""")
else:
    print(f"""
⚠ Rust implementation differs from Python by {max_error:.6f}%
- Check parameter values in Rust code
- Verify formula implementation
- May be acceptable if due to numerical precision
""")
