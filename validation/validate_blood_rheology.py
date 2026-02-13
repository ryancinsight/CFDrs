#!/usr/bin/env python3
"""
Blood Rheology Validation

Validates pycfdrs blood rheology models (Casson, Cross, Carreau-Yasuda)
against expected shear-thinning behavior and literature parameters.
"""

import numpy as np
import pycfdrs
import matplotlib.pyplot as plt
from pathlib import Path

def validate_casson_blood():
    """Validate Casson blood model shear-thinning behavior."""
    print(f"\n{'='*70}")
    print(f"Casson Blood Model Validation")
    print(f"{'='*70}\n")
    
    blood = pycfdrs.CassonBlood()
    
    # Shear rate range (typical for blood flow: 1-1000 s⁻¹)
    gamma_dot = np.logspace(0, 3, 100)  # 1 to 1000 s⁻¹
    
    # Compute viscosity at each shear rate
    mu = np.array([blood.apparent_viscosity(g) for g in gamma_dot])
    
    # Expected behavior
    mu_inf = 3.5e-3  # High shear viscosity (Pa·s)
    
    print(f"Casson Parameters:")
    print(f"  τ_y = {blood.yield_stress() * 1000:.2f} mPa")
    print(f"  μ_∞ = {blood.viscosity_high_shear() * 1000:.2f} mPa·s\n")
    
    print(f"Shear-Rate Dependent Viscosity:")
    print(f"  μ(1 s⁻¹)    = {mu[0] * 1000:.2f} mPa·s")
    print(f"  μ(10 s⁻¹)   = {mu[10] * 1000:.2f} mPa·s")
    print(f"  μ(100 s⁻¹)  = {mu[30] * 1000:.2f} mPa·s")
    print(f"  μ(1000 s⁻¹) = {mu[-1] * 1000:.2f} mPa·s\n")
    
    # Validate shear-thinning
    is_shear_thinning = np.all(np.diff(mu) <= 0.0001)  # Monotonic decrease (with small tolerance)
    converges_to_inf = abs(mu[-1] - mu_inf) / mu_inf < 0.10  # Within 10% of μ_∞
    
    print(f"Validation Checks:")
    print(f"  Shear-thinning:      {'[OK]' if is_shear_thinning else '[FAIL]'}")
    print(f"  Converges to mu_inf:    {'[OK]' if converges_to_inf else '[FAIL]'}")
    
    passed = is_shear_thinning and converges_to_inf
    print(f"\nStatus: {'PASSED' if passed else 'FAILED'}\n")
    
    # Plot viscosity vs shear rate
    plt.figure(figsize=(10, 6))
    plt.loglog(gamma_dot, mu * 1000, 'b-', linewidth=2, label='Casson Model')
    plt.axhline(y=mu_inf * 1000, color='r', linestyle='--', label=f'μ_∞ = {mu_inf * 1000:.1f} mPa·s')
    plt.xlabel('Shear Rate γ̇ (s⁻¹)', fontsize=12)
    plt.ylabel('Apparent Viscosity μ (mPa·s)', fontsize=12)
    plt.title('Casson Blood: Shear-Thinning Behavior', fontsize=14, fontweight='bold')
    plt.legend(fontsize=11)
    plt.grid(True, which='both', alpha=0.3)
    plt.tight_layout()
    
    output_dir = Path(__file__).parent / "results"
    output_dir.mkdir(exist_ok=True)
    plt.savefig(output_dir / "casson_rheology.png", dpi=150)
    print(f"Plot saved: {output_dir / 'casson_rheology.png'}\n")
    
    return passed

def validate_cross_blood():
    """Validate Cross blood model."""
    print(f"\n{'='*70}")
    print(f"Cross Blood Model Validation")
    print(f"{'='*70}\n")
    
    blood = pycfdrs.CrossBlood()
    
    gamma_dot = np.logspace(-1, 3, 100)
    mu = np.array([blood.apparent_viscosity(g) for g in gamma_dot])
    
    # Cross model parameters (typical for blood)
    mu_0 = 56e-3  # Zero-shear viscosity (Pa·s)
    mu_inf = 3.5e-3  # Infinite-shear viscosity (Pa·s)
    
    print(f"Cross Model Limits:")
    print(f"  μ_0:  {mu[0] * 1000:.2f} mPa·s (expected ≈ {mu_0 * 1000:.0f})")
    print(f"  μ_∞:  {mu[-1] * 1000:.2f} mPa·s (expected ≈ {mu_inf * 1000:.1f})\n")
    
    # Validate limits
    mu_0_match = abs(mu[0] - mu_0) / mu_0 < 0.20
    mu_inf_match = abs(mu[-1] - mu_inf) / mu_inf < 0.20
    
    print(f"Validation:")
    print(f"  mu_0 match:  {'[OK]' if mu_0_match else '[FAIL]'}")
    print(f"  mu_inf match:  {'[OK]' if mu_inf_match else '[FAIL]'}")
    
    passed = mu_0_match and mu_inf_match
    print(f"\nStatus: {'PASSED' if passed else 'FAILED'}\n")
    
    return passed

def compare_models():
    """Compare all blood models."""
    print(f"\n{'='*70}")
    print(f"Comparing Blood Rheology Models")
    print(f"{'='*70}\n")
    
    gamma_dot = np.logspace(0, 3, 100)
    
    casson = pycfdrs.CassonBlood()
    cross = pycfdrs.CrossBlood()
    carreau = pycfdrs.CarreauYasudaBlood()
    
    mu_casson = np.array([casson.apparent_viscosity(g) for g in gamma_dot])
    mu_cross = np.array([cross.apparent_viscosity(g) for g in gamma_dot])
    mu_carreau = np.array([carreau.apparent_viscosity(g) for g in gamma_dot])
    
    plt.figure(figsize=(12, 7))
    plt.loglog(gamma_dot, mu_casson * 1000, 'b-', linewidth=2, label='Casson')
    plt.loglog(gamma_dot, mu_cross * 1000, 'r--', linewidth=2, label='Cross')
    plt.loglog(gamma_dot, mu_carreau * 1000, 'g-.', linewidth=2, label='Carreau-Yasuda')
    plt.xlabel('Shear Rate γ̇ (s⁻¹)', fontsize=12)
    plt.ylabel('Apparent Viscosity μ (mPa·s)', fontsize=12)
    plt.title('Blood Rheology Models Comparison', fontsize=14, fontweight='bold')
    plt.legend(fontsize=11)
    plt.grid(True, which='both', alpha=0.3)
    plt.tight_layout()
    
    output_dir = Path(__file__).parent / "results"
    output_dir.mkdir(exist_ok=True)
    plt.savefig(output_dir / "blood_models_comparison.png", dpi=150)
    print(f"Comparison plot saved: {output_dir / 'blood_models_comparison.png'}\n")

if __name__ == "__main__":
    casson_pass = validate_casson_blood()
    cross_pass = validate_cross_blood()
    compare_models()
    
    print(f"\n{'='*70}")
    print(f"Summary:")
    print(f"  Casson:  {'PASSED' if casson_pass else 'FAILED'}")
    print(f"  Cross:   {'PASSED' if cross_pass else 'FAILED'}")
    print(f"{'='*70}\n")
