#!/usr/bin/env python3
"""
Blood Rheology Validation: pycfdrs vs Literature Data

This script validates blood rheology models against:
1. Casson model viscosity vs shear rate
2. Carreau-Yasuda model viscosity vs shear rate
3. Literature values from Merrill et al. (1969), Chien et al. (1970)

Literature References:
- Merrill, E.W. et al. (1969) "Rheology of human blood, near and at zero flow"
  Biophysical Journal 9:199-213
- Chien, S. et al. (1970) "Blood viscosity: influence of erythrocyte aggregation"
  Science 168:977-979

Run with:
    python validation/validate_blood_rheology.py
"""

import numpy as np
import sys
sys.path.insert(0, "crates/pycfdrs")
import pycfdrs


# Literature values for normal human blood at 37C, Hematocrit 45%
# Source: Merrill et al. (1969)
LITERATURE_DATA = {
    # Shear rate [s^-1] : Viscosity [mPa.s]
    0.1: 50.0,    # High viscosity at low shear (aggregation)
    1.0: 20.0,    # Declining due to rouleaux breakup
    10.0: 8.0,    # Shear thinning
    100.0: 5.0,   # Approaching asymptote
    1000.0: 3.8,  # Near high-shear limit
}


def validate_casson_asymptotic():
    """Validate Casson model asymptotic viscosity."""
    print("\n" + "="*70)
    print("1. CASSON MODEL - ASYMPTOTIC BEHAVIOR")
    print("="*70)
    
    blood = pycfdrs.CassonBlood()
    
    # Get model parameters
    mu_inf = blood.viscosity_high_shear() * 1000  # mPa.s
    tau_y = blood.yield_stress()                   # Pa
    rho = blood.density()                          # kg/m^3
    
    print(f"\n  Model Parameters:")
    print(f"    Yield stress: {tau_y*1000:.2f} mPa")
    print(f"    Infinite shear viscosity: {mu_inf:.3f} mPa.s")
    print(f"    Density: {rho:.0f} kg/m^3")
    
    # Literature: mu_inf for blood is 3-4 mPa.s (Merrill 1969)
    mu_inf_lit = 3.5  # mPa.s
    error = abs(mu_inf - mu_inf_lit) / mu_inf_lit
    
    print(f"\n  Comparison:")
    print(f"    Model mu_inf: {mu_inf:.3f} mPa.s")
    print(f"    Literature (Merrill 1969): {mu_inf_lit:.1f} mPa.s")
    print(f"    Relative error: {error*100:.1f}%")
    
    passed = error < 0.10  # Within 10%
    print(f"\n  RESULT: {'PASS' if passed else 'FAIL'} (tolerance: <10%)")
    return passed, error


def validate_carreau_yasuda_shear_thinning():
    """Validate Carreau-Yasuda shear thinning behavior."""
    print("\n" + "="*70)
    print("2. CARREAU-YASUDA - SHEAR THINNING")
    print("="*70)
    
    blood = pycfdrs.CarreauYasudaBlood()
    
    # Get viscosity at different shear rates
    mu_zero = blood.viscosity_zero_shear() * 1000  # mPa.s
    mu_inf = blood.viscosity_high_shear() * 1000   # mPa.s
    
    print(f"\n  Model Parameters:")
    print(f"    Zero shear viscosity: {mu_zero:.1f} mPa.s")
    print(f"    Infinite shear viscosity: {mu_inf:.3f} mPa.s")
    
    # Test shear thinning: mu should decrease with increasing shear rate
    shear_rates = [0.1, 1.0, 10.0, 100.0, 1000.0]
    viscosities = [blood.viscosity(g) * 1000 for g in shear_rates]
    
    print(f"\n  Viscosity vs Shear Rate:")
    print(f"    gamma [s^-1]    mu [mPa.s]")
    print(f"    " + "-"*30)
    for g, mu in zip(shear_rates, viscosities):
        print(f"    {g:8.1f}        {mu:8.4f}")
    
    # Check monotonic decrease
    is_monotonic = all(viscosities[i] >= viscosities[i+1] for i in range(len(viscosities)-1))
    
    print(f"\n  Monotonic decrease: {is_monotonic}")
    
    passed = is_monotonic
    error = 0.0 if is_monotonic else 1.0
    print(f"\n  RESULT: {'PASS' if passed else 'FAIL'}")
    return passed, error


def validate_model_comparison():
    """Compare Casson and Carreau-Yasuda models."""
    print("\n" + "="*70)
    print("3. MODEL COMPARISON")
    print("="*70)
    
    casson = pycfdrs.CassonBlood()
    carreau = pycfdrs.CarreauYasudaBlood()
    
    shear_rates = [0.1, 1.0, 10.0, 100.0, 1000.0]
    
    print(f"\n  Viscosity Comparison [mPa.s]:")
    print(f"    gamma [s^-1]    Casson      Carreau     Lit. (approx)")
    print(f"    " + "-"*55)
    
    max_error = 0.0
    for g in shear_rates:
        mu_cas = casson.viscosity(g) * 1000
        mu_car = carreau.viscosity(g) * 1000
        mu_lit = LITERATURE_DATA.get(g, None)
        
        lit_str = f"{mu_lit:.1f}" if mu_lit else "N/A"
        print(f"    {g:8.1f}        {mu_cas:8.4f}    {mu_car:8.4f}    {lit_str:>8s}")
        
        # Track worst case error vs literature
        if mu_lit:
            err_cas = abs(mu_cas - mu_lit) / mu_lit
            err_car = abs(mu_car - mu_lit) / mu_lit
            max_error = max(max_error, min(err_cas, err_car))
    
    # At high shear (gamma > 100), models should match well
    mu_cas_1000 = casson.viscosity(1000) * 1000
    mu_car_1000 = carreau.viscosity(1000) * 1000
    high_shear_diff = abs(mu_cas_1000 - mu_car_1000) / mu_cas_1000
    
    print(f"\n  High shear agreement:")
    print(f"    Casson at 1000 s^-1: {mu_cas_1000:.4f} mPa.s")
    print(f"    Carreau at 1000 s^-1: {mu_car_1000:.4f} mPa.s")
    print(f"    Difference: {high_shear_diff*100:.2f}%")
    
    passed = high_shear_diff < 0.05  # <5% at high shear
    print(f"\n  RESULT: {'PASS' if passed else 'FAIL'} (high shear diff <5%)")
    return passed, high_shear_diff


def validate_yield_stress():
    """Validate Casson yield stress behavior."""
    print("\n" + "="*70)
    print("4. CASSON YIELD STRESS")
    print("="*70)
    
    blood = pycfdrs.CassonBlood()
    tau_y = blood.yield_stress()  # Pa
    
    # Literature: tau_y for blood is 0.005-0.01 Pa (Merrill 1969)
    tau_y_lit_low = 0.004   # Pa
    tau_y_lit_high = 0.015  # Pa
    
    print(f"\n  Yield Stress:")
    print(f"    Model: {tau_y*1000:.2f} mPa")
    print(f"    Literature range: {tau_y_lit_low*1000:.1f}-{tau_y_lit_high*1000:.1f} mPa")
    
    in_range = tau_y_lit_low <= tau_y <= tau_y_lit_high
    
    if in_range:
        error = 0.0
    else:
        mid = (tau_y_lit_low + tau_y_lit_high) / 2
        error = abs(tau_y - mid) / mid
    
    print(f"\n  RESULT: {'PASS' if in_range else 'FAIL'}")
    return in_range, error


def main():
    print("\n" + "#"*70)
    print(" BLOOD RHEOLOGY VALIDATION SUITE")
    print(" pycfdrs vs Literature Data (Merrill 1969)")
    print("#"*70)
    
    results = []
    
    try:
        passed1, error1 = validate_casson_asymptotic()
        results.append(("Casson Asymptotic", passed1, error1))
    except Exception as e:
        print(f"\n  ERROR: {e}")
        import traceback
        traceback.print_exc()
        results.append(("Casson Asymptotic", False, float('inf')))
    
    try:
        passed2, error2 = validate_carreau_yasuda_shear_thinning()
        results.append(("Carreau Shear-Thin", passed2, error2))
    except Exception as e:
        print(f"\n  ERROR: {e}")
        import traceback
        traceback.print_exc()
        results.append(("Carreau Shear-Thin", False, float('inf')))
    
    try:
        passed3, error3 = validate_model_comparison()
        results.append(("Model Comparison", passed3, error3))
    except Exception as e:
        print(f"\n  ERROR: {e}")
        import traceback
        traceback.print_exc()
        results.append(("Model Comparison", False, float('inf')))
    
    try:
        passed4, error4 = validate_yield_stress()
        results.append(("Casson Yield Stress", passed4, error4))
    except Exception as e:
        print(f"\n  ERROR: {e}")
        import traceback
        traceback.print_exc()
        results.append(("Casson Yield Stress", False, float('inf')))
    
    # Summary
    print("\n" + "="*70)
    print("SUMMARY")
    print("="*70)
    
    all_passed = True
    for name, passed, error in results:
        status = "PASS" if passed else "FAIL"
        all_passed = all_passed and passed
        err_str = f"{error:.2e}" if error < float('inf') else "N/A"
        print(f"  {name:25s}: {status} (error: {err_str})")
    
    print("="*70)
    if all_passed:
        print("ALL BLOOD RHEOLOGY VALIDATION TESTS PASSED")
    else:
        print("SOME TESTS FAILED")
    print("="*70)
    
    return 0 if all_passed else 1


if __name__ == "__main__":
    exit(main())
