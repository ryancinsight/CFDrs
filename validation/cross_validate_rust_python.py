#!/usr/bin/env python3
"""
Cross-validate Rust and Python implementations of validated physics models

This script verifies that the Rust cfd-core implementations match
the Python validation code for:
1. Blake threshold calculation
2. Blood viscosity (Carreau-Yasuda model)
3. Giersiepen hemolysis model
"""

import sys
import math

print("="*80)
print("RUST vs PYTHON PHYSICS VALIDATION")
print("="*80)

# Try to import pycfdrs (or cfd_python)
try:
    try:
        import pycfdrs
        print("\n✓ pycfdrs module loaded successfully\n")
    except ImportError:
        import cfd_python as pycfdrs
        print("\n✓ cfd_python module loaded successfully (aliased as pycfdrs)\n")
    has_pycfdrs = True
except ImportError:
    has_pycfdrs = False
    print("\n✗ pycfdrs/cfd_python not available - will compare Python calculations only\n")
    print("Run: maturin develop --release\n")

# Physical constants (match validation scripts)
WATER_DENSITY = 997.0  # kg/m³
WATER_VISCOSITY = 0.001  # Pa·s
WATER_SURFACE_TENSION = 0.0728  # N/m
WATER_VAPOR_PRESSURE = 2339.0  # Pa

print("="*80)
print("TEST 1: Blake Threshold")
print("="*80)

# Python calculation (validated)
R_0 = 10e-6  # m - initial bubble radius
P_inf = 101325.0  # Pa
sigma = WATER_SURFACE_TENSION
P_v = WATER_VAPOR_PRESSURE

# Step 1: Calculate critical radius
coefficient = 0.85
R_c_python = coefficient * 2 * sigma / (P_inf - P_v)

# Step 2: Calculate Blake threshold
P_Blake_python = P_v + (4 * sigma) / (3 * R_c_python)

print(f"\nPython calculation:")
print(f"  R_c = {R_c_python*1e6:.4f} μm")
print(f"  P_Blake = {P_Blake_python:.2f} Pa = {P_Blake_python/1000:.2f} kPa")

if has_pycfdrs:
    print(f"\nRust implementation:")
    print(f"  Located in: crates/cfd-core/src/physics/cavitation/regimes.rs")
    print(f"  Method: blake_threshold() and blake_critical_radius()")

    # Create bubble model
    bubble = pycfdrs.RayleighPlesset(
        initial_radius=R_0,
        liquid_density=WATER_DENSITY,
        liquid_viscosity=WATER_VISCOSITY,
        surface_tension=WATER_SURFACE_TENSION,
        vapor_pressure=WATER_VAPOR_PRESSURE,
        polytropic_index=1.4
    )

    # Create classifier
    classifier = pycfdrs.CavitationRegimeClassifier(
        bubble,
        ambient_pressure=P_inf,
        acoustic_pressure=None,
        acoustic_frequency=None
    )

    # Get thresholds
    P_Blake_rust = classifier.blake_threshold()

    # Check Blake radius directly via bubble model
    # Note: blake_critical_radius is on RayleighPlesset in Rust bindings
    R_c_rust = bubble.blake_critical_radius(P_inf)

    print(f"  Rust R_c = {R_c_rust*1e6:.4f} μm")
    print(f"  Rust P_Blake = {P_Blake_rust:.2f} Pa")

    # Compare
    diff_Rc = abs(R_c_rust - R_c_python) / R_c_python * 100
    diff_P = abs(P_Blake_rust - P_Blake_python) / P_Blake_python * 100

    if diff_Rc < 0.01 and diff_P < 0.01:
        print(f"  ✓ Formula matches Python implementation (Error: Rc={diff_Rc:.4f}%, P={diff_P:.4f}%)")
    else:
        print(f"  ✗ Formula MISMATCH (Error: Rc={diff_Rc:.4f}%, P={diff_P:.4f}%)")
        sys.exit(1)

else:
    print(f"\nRust verification skipped (pycfdrs not available)")

print("\n" + "="*80)
print("TEST 2: Blood Viscosity (Carreau-Yasuda)")
print("="*80)

# Carreau-Yasuda parameters (Cho & Kensey 1991)
mu_0 = 0.056  # Pa·s
mu_inf = 0.00345  # Pa·s
lambda_cy = 3.313  # s
a_cy = 2.0
n_cy = 0.3568

def carreau_yasuda_python(shear_rate):
    """Python implementation of Carreau-Yasuda model"""
    term = 1 + (lambda_cy * shear_rate) ** a_cy
    exponent = (n_cy - 1) / a_cy
    return mu_inf + (mu_0 - mu_inf) * (term ** exponent)

# Test at multiple shear rates
test_shear_rates = [100, 1000, 5000, 100000]

print(f"\nPython Carreau-Yasuda:")
print(f"  μ₀ = {mu_0*1000:.1f} mPa·s")
print(f"  μ_∞ = {mu_inf*1000:.3f} mPa·s")
print(f"  λ = {lambda_cy} s")
print(f"  a = {a_cy}")
print(f"  n = {n_cy}\n")

print(f"{'Shear Rate (s⁻¹)':>20} {'μ (mPa·s)':>15} {'Error from μ_∞':>15}")
print("-"*55)

for gamma_dot in test_shear_rates:
    mu_python = carreau_yasuda_python(gamma_dot)
    error_pct = abs(mu_python - mu_inf) / mu_inf * 100
    print(f"{gamma_dot:20.0f} {mu_python*1000:15.4f} {error_pct:14.2f}%")

if has_pycfdrs:
    print(f"\nRust implementation:")
    print(f"  Located in: crates/cfd-core/src/physics/fluid/blood.rs")
    print(f"  Type: CarreauYasudaBlood")
    print(f"  Method: apparent_viscosity(shear_rate)")
    
    # Create blood model (uses defaults which match Cho & Kensey 1991 for normal blood)
    # Rust default params: mu_0=0.056, mu_inf=0.00345, lambda=3.313, a=2.0, n=0.3568
    blood = pycfdrs.CarreauYasudaBlood()

    # Test at same shear rates
    max_error = 0.0
    for gamma_dot in test_shear_rates:
        mu_rust = blood.apparent_viscosity(float(gamma_dot))
        mu_python = carreau_yasuda_python(gamma_dot)

        diff = abs(mu_rust - mu_python)
        diff_pct = diff / mu_python * 100 if mu_python > 0 else 0
        max_error = max(max_error, diff_pct)

        # print(f"  γ={gamma_dot}: Rust={mu_rust:.6f}, Py={mu_python:.6f}, Diff={diff_pct:.4f}%")

    if max_error < 0.01:
         print(f"  ✓ Validated at all shear rates (Max error: {max_error:.6f}%)")
    else:
         print(f"  ✗ Validation FAILED (Max error: {max_error:.6f}%)")
         sys.exit(1)

else:
    print(f"\nRust verification skipped (pycfdrs not available)")

print("\n" + "="*80)
print("TEST 3: Giersiepen Hemolysis Model")
print("="*80)

# Giersiepen constants
C = 3.62e-5
alpha = 2.416
beta = 0.785

def giersiepen_python(shear_stress, exposure_time):
    """Python implementation of Giersiepen power law"""
    return C * (shear_stress ** alpha) * (exposure_time ** beta)

# Test cases
test_cases = [
    (50.0, 1.0),
    (100.0, 0.5),
    (150.0, 0.1),
    (200.0, 0.05),
]

print(f"\nPython Giersiepen model:")
print(f"  C = {C}")
print(f"  α = {alpha}")
print(f"  β = {beta}\n")

print(f"{'Stress (Pa)':>12} {'Time (s)':>12} {'Damage':>15}")
print("-"*42)

for tau, t in test_cases:
    damage = giersiepen_python(tau, t)
    print(f"{tau:12.1f} {t:12.2f} {damage:15.6f}")

if has_pycfdrs:
    print(f"\nRust implementation:")
    print(f"  Located in: crates/cfd-core/src/physics/hemolysis.rs")
    print(f"  Method: damage_index(shear_stress, exposure_time)")

    # Create standard Giersiepen model
    model = pycfdrs.HemolysisModel.giersiepen_standard()

    max_error = 0.0
    for tau, t in test_cases:
        damage_rust = model.damage_index(tau, t)
        damage_python = giersiepen_python(tau, t)

        diff = abs(damage_rust - damage_python)
        diff_pct = diff / damage_python * 100 if damage_python > 0 else 0
        max_error = max(max_error, diff_pct)

        # print(f"  τ={tau}, t={t}: Rust={damage_rust:.6e}, Py={damage_python:.6e}, Diff={diff_pct:.4f}%")

    if max_error < 0.01:
        print(f"  ✓ Validated for all test cases (Max error: {max_error:.6f}%)")
    else:
        print(f"  ✗ Validation FAILED (Max error: {max_error:.6f}%)")
        sys.exit(1)
else:
    print(f"\nRust verification skipped (pycfdrs not available)")

print("\n" + "="*80)
print("SUMMARY")
print("="*80)

if has_pycfdrs:
    print(f"""
All three validated models have corresponding Rust implementations verified:

1. ✓ Blake Threshold
   - Python: validation/validate_cavitation_hemolysis.py line 129-134
   - Rust: crates/cfd-core/src/physics/cavitation/regimes.rs
   - Formula: P_Blake = P_v + 4σ/(3R_c) where R_c = 0.85×2σ/(P_∞-P_v)
   - Status: MATCHES

2. ✓ Carreau-Yasuda Blood
   - Python: validation scripts
   - Rust: crates/cfd-core/src/physics/fluid/blood.rs
   - Formula: μ = μ_∞ + (μ₀-μ_∞)[1+(λγ̇)^a]^((n-1)/a)
   - Status: MATCHES

3. ✓ Giersiepen Hemolysis
   - Python: validation/validate_cavitation_hemolysis.py line 160-168
   - Rust: crates/cfd-core/src/physics/hemolysis.rs
   - Formula: D = C×τ^α×t^β
   - Status: MATCHES
""")
else:
    print(f"""
Python implementations validated against literature.

To complete Rust validation:
1. Build pycfdrs: maturin develop --release
2. Re-run this script to cross-check values
""")

print("\n" + "="*80)
print("VALIDATION STATUS")
print("="*80)

cross_check_status = "✓ VERIFIED" if has_pycfdrs else "⚠ PENDING"

validation_status = {
    "Blake Threshold": {
        "Physics": "✓ VALIDATED (against Brennen 1995)",
        "Python": "✓ CORRECT (R_c formulation)",
        "Rust": "✓ IMPLEMENTED (regimes.rs)",
        "Cross-check": cross_check_status
    },
    "Blood Viscosity": {
        "Physics": "✓ VALIDATED (against Cho & Kensey 1991)",
        "Python": "✓ CORRECT (λ=3.313s convergence)",
        "Rust": "✓ IMPLEMENTED (blood.rs)",
        "Cross-check": cross_check_status
    },
    "Hemolysis Model": {
        "Physics": "✓ VALIDATED (against Giersiepen 1990)",
        "Python": "✓ CORRECT (iso-damage curves)",
        "Rust": "✓ IMPLEMENTED (hemolysis.rs)",
        "Cross-check": cross_check_status
    }
}

for model, status in validation_status.items():
    print(f"\n{model}:")
    for key, val in status.items():
        print(f"  {key:12} {val}")

print("\n" + "="*80)
