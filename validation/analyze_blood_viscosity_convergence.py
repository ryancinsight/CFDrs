#!/usr/bin/env python3
"""
Deep investigation of Carreau-Yasuda blood viscosity convergence
"""

import math

# Cho & Kensey (1991) parameters
mu_0 = 0.056  # Pa·s - zero shear viscosity
mu_inf = 0.00345  # Pa·s - infinite shear viscosity
lambda_val = 3.313  # s - relaxation time
n = 0.3568  # power-law index
a = 2.0  # transition parameter

print("="*80)
print("CARREAU-YASUDA BLOOD VISCOSITY CONVERGENCE ANALYSIS")
print("="*80)
print(f"\nModel parameters (Cho & Kensey, 1991):")
print(f"  μ₀ = {mu_0*1000:.3f} mPa·s (zero shear)")
print(f"  μ_∞ = {mu_inf*1000:.3f} mPa·s (infinite shear)")
print(f"  λ = {lambda_val} s")
print(f"  n = {n}")
print(f"  a = {a}")

print(f"\n{'γ̇ (s⁻¹)':>12} {'μ (mPa·s)':>12} {'Error from μ_∞':>18} {'Converged?':>12}")
print("-"*60)

shear_rates = [0.1, 1, 10, 100, 1000, 5000, 10000, 50000, 100000, 1000000]

for gamma in shear_rates:
    # Carreau-Yasuda model: μ = μ_∞ + (μ_0 - μ_∞) · [1 + (λγ̇)^a]^((n-1)/a)
    lambda_gamma = lambda_val * gamma
    lambda_gamma_a = lambda_gamma ** a
    bracketed = 1.0 + lambda_gamma_a
    exponent = (n - 1.0) / a
    shear_factor = bracketed ** exponent
    
    mu = mu_inf + (mu_0 - mu_inf) * shear_factor
    error_pct = abs(mu - mu_inf) / mu_inf * 100
    
    converged_5pct = "YES ✓" if error_pct < 5.0 else "NO"
    converged_10pct = "YES ✓" if error_pct < 10.0 else "NO"
    
    print(f"{gamma:12.1f} {mu*1000:12.4f} {error_pct:15.2f}% {converged_5pct:>12}")

print("\n" + "="*80)
print("ANALYSIS")
print("="*80)

# Find shear rate for 5% convergence
target_error = 0.05
# μ = μ_∞ + (μ_0 - μ_∞) · [1 + (λγ̇)^a]^((n-1)/a)
# For 5% error: μ = 1.05 * μ_∞
# 1.05 * μ_∞ = μ_∞ + (μ_0 - μ_∞) · [1 + (λγ̇)^a]^((n-1)/a)
# 0.05 * μ_∞ / (μ_0 - μ_∞) = [1 + (λγ̇)^a]^((n-1)/a)

shear_factor_target = 0.05 * mu_inf / (mu_0 - mu_inf)
# [1 + (λγ̇)^a]^((n-1)/a) = shear_factor_target
# 1 + (λγ̇)^a = shear_factor_target^(a/(n-1))
# (λγ̇)^a = shear_factor_target^(a/(n-1)) - 1
# λγ̇ = [shear_factor_target^(a/(n-1)) - 1]^(1/a)
# γ̇ = [shear_factor_target^(a/(n-1)) - 1]^(1/a) / λ

exponent_inv = a / (n - 1.0)
bracketed_target = shear_factor_target ** exponent_inv
lambda_gamma_a_target = bracketed_target - 1.0

if lambda_gamma_a_target > 0:
    lambda_gamma_target = lambda_gamma_a_target ** (1.0 / a)
    gamma_target_5pct = lambda_gamma_target / lambda_val
    print(f"\nFor 5% convergence to μ_∞:")
    print(f"  Required shear rate: γ̇ ≈ {gamma_target_5pct:.0f} s⁻¹")
else:
    print(f"\nFor 5% convergence to μ_∞:")
    print(f"  Model NEVER converges within 5% for these parameters!")

# Check if λ value is correct from literature
print(f"\n{'='*80}")
print("LITERATURE REVIEW - Cho & Kensey (1991)")
print("="*80)
print("""
The Carreau-Yasuda parameters from Cho & Kensey (1991):
- μ₀ = 0.056 Pa·s
- μ_∞ = 0.00345 Pa·s  
- λ = 3.313 s  ← This is VERY LARGE!
- n = 0.3568
- a = 2.0

The relaxation time λ = 3.313 s means the model predicts:
- Even at γ̇ = 1000 s⁻¹, we have λγ̇ = 3313
- This gives (λγ̇)^a = 3313^2 ≈ 1.1×10⁷
- The model is still FAR from infinite-shear regime

CONCLUSION: The 5% convergence criterion is UNREALISTIC for these
literature-validated parameters. The Carreau-Yasuda model with 
λ = 3.313 s predicts VERY SLOW convergence to μ_∞.

At physiologically relevant shear rates (10-1000 s⁻¹), the model
is designed to be in the transitional regime, NOT the infinite-shear
regime.
""")

print("="*80)
print("RECOMMENDATIONS")
print("="*80)
print("""
1. ACCEPT that 8-10% error at γ̇ = 1000 s⁻¹ is CORRECT MODEL BEHAVIOR
   - This is what the Cho & Kensey parameters predict
   - Venturi throats (~1000 s⁻¹) are still in transitional regime
   
2. TEST at γ̇ = 100,000 s⁻¹ for <5% convergence
   - This is beyond physiological range
   - Only for validating infinite-shear limit
   
3. OR test convergence at realistic shear rates with appropriate tolerance:
   - γ̇ = 1000 s⁻¹ → expect 8-10% from μ_∞ (CORRECT)
   - γ̇ = 100 s⁻¹ → expect 25-35% from μ_∞ (CORRECT)
   
The original test tolerance of 5% was TOO STRICT for the 
literature-validated Carreau-Yasuda parameters at physiological shear rates.
""")
