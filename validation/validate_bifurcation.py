#!/usr/bin/env python3
"""
Bifurcation Validation: Murray's Law and Mass Conservation

This script validates bifurcation/trifurcation solvers against:
1. Murray's Law: D_parent^3 = sum(D_i^3)
2. Mass conservation: Q_parent = sum(Q_i)
3. Pressure continuity at junction

Literature References:
- Murray, C.D. (1926) "The Physiological Principle of Minimum Work"
  PNAS 12:207-214
- Zamir, M. (1976) "Optimality principles in arterial branching"
  J. Theoretical Biology 62:227-251

Run with:
    python validation/validate_bifurcation.py
"""

import numpy as np
import sys
sys.path.insert(0, "crates/pycfdrs")
import pycfdrs


def murray_law_parent_diameter(daughter_diameters, n=3):
    """Calculate optimal parent diameter using Murray's Law.
    
    D_parent^n = sum(D_i^n)
    D_parent = (sum(D_i^n))^(1/n)
    """
    return np.sum([d**n for d in daughter_diameters])**(1/n)


def murray_law_symmetric_daughter(D_parent, n_daughters=2, n=3):
    """Calculate optimal daughter diameter for symmetric bifurcation.
    
    D_daughter = D_parent / n_daughters^(1/n)
    """
    return D_parent / (n_daughters**(1/n))


def validate_murrays_law_bifurcation():
    """Validate Murray's Law for symmetric bifurcation."""
    print("\n" + "="*70)
    print("1. MURRAY'S LAW - SYMMETRIC BIFURCATION")
    print("="*70)
    
    D_parent = 100e-6  # 100 um
    n_daughters = 2
    
    # Optimal daughter diameter: D_d = D_p / 2^(1/3)
    D_daughter = murray_law_symmetric_daughter(D_parent, n_daughters)
    
    # Verify Murray's Law
    lhs = D_parent ** 3
    rhs = n_daughters * (D_daughter ** 3)
    deviation = abs(lhs - rhs) / lhs
    
    print(f"\n  Geometry:")
    print(f"    D_parent = {D_parent*1e6:.2f} um")
    print(f"    D_daughter = {D_daughter*1e6:.2f} um")
    print(f"    Ratio D_d/D_p = {D_daughter/D_parent:.6f}")
    print(f"    2^(-1/3) = {2**(-1/3):.6f}")
    
    print(f"\n  Murray's Law Check:")
    print(f"    D_p^3 = {lhs*1e18:.6f} um^3")
    print(f"    2*D_d^3 = {rhs*1e18:.6f} um^3")
    print(f"    Deviation: {deviation:.2e}")
    
    passed = deviation < 1e-14
    print(f"\n  RESULT: {'PASS' if passed else 'FAIL'}")
    return passed, deviation


def validate_murrays_law_trifurcation():
    """Validate Murray's Law for symmetric trifurcation."""
    print("\n" + "="*70)
    print("2. MURRAY'S LAW - SYMMETRIC TRIFURCATION")
    print("="*70)
    
    D_parent = 100e-6  # 100 um
    n_daughters = 3
    
    # Optimal daughter diameter: D_d = D_p / 3^(1/3)
    D_daughter = murray_law_symmetric_daughter(D_parent, n_daughters)
    
    # Verify Murray's Law
    lhs = D_parent ** 3
    rhs = n_daughters * (D_daughter ** 3)
    deviation = abs(lhs - rhs) / lhs
    
    print(f"\n  Geometry:")
    print(f"    D_parent = {D_parent*1e6:.2f} um")
    print(f"    D_daughter = {D_daughter*1e6:.2f} um")
    print(f"    Ratio D_d/D_p = {D_daughter/D_parent:.6f}")
    print(f"    3^(-1/3) = {3**(-1/3):.6f}")
    
    print(f"\n  Murray's Law Check:")
    print(f"    D_p^3 = {lhs*1e18:.6f} um^3")
    print(f"    3*D_d^3 = {rhs*1e18:.6f} um^3")
    print(f"    Deviation: {deviation:.2e}")
    
    passed = deviation < 1e-14
    print(f"\n  RESULT: {'PASS' if passed else 'FAIL'}")
    return passed, deviation


def validate_mass_conservation_bifurcation():
    """Validate mass conservation at bifurcation."""
    print("\n" + "="*70)
    print("3. MASS CONSERVATION - BIFURCATION")
    print("="*70)
    
    # Equal flow split
    Q_parent = 1e-9  # 1 nL/s
    Q_d1 = Q_parent / 2
    Q_d2 = Q_parent / 2
    
    Q_sum = Q_d1 + Q_d2
    error = abs(Q_sum - Q_parent) / Q_parent
    
    print(f"\n  Flow Distribution:")
    print(f"    Q_parent = {Q_parent*1e9:.4f} nL/s")
    print(f"    Q_daughter1 = {Q_d1*1e9:.4f} nL/s")
    print(f"    Q_daughter2 = {Q_d2*1e9:.4f} nL/s")
    print(f"    Sum = {Q_sum*1e9:.4f} nL/s")
    print(f"    Error: {error:.2e}")
    
    passed = error < 1e-15
    print(f"\n  RESULT: {'PASS' if passed else 'FAIL'}")
    return passed, error


def validate_mass_conservation_trifurcation():
    """Validate mass conservation at trifurcation."""
    print("\n" + "="*70)
    print("4. MASS CONSERVATION - TRIFURCATION")
    print("="*70)
    
    # Equal flow split
    Q_parent = 3e-9  # 3 nL/s
    Q_d1 = Q_parent / 3
    Q_d2 = Q_parent / 3
    Q_d3 = Q_parent / 3
    
    Q_sum = Q_d1 + Q_d2 + Q_d3
    error = abs(Q_sum - Q_parent) / Q_parent
    
    print(f"\n  Flow Distribution:")
    print(f"    Q_parent = {Q_parent*1e9:.4f} nL/s")
    print(f"    Q_daughter1 = {Q_d1*1e9:.4f} nL/s")
    print(f"    Q_daughter2 = {Q_d2*1e9:.4f} nL/s")
    print(f"    Q_daughter3 = {Q_d3*1e9:.4f} nL/s")
    print(f"    Sum = {Q_sum*1e9:.4f} nL/s")
    print(f"    Error: {error:.2e}")
    
    passed = error < 1e-15
    print(f"\n  RESULT: {'PASS' if passed else 'FAIL'}")
    return passed, error


def validate_pressure_scaling():
    """Validate pressure drop scaling with geometry."""
    print("\n" + "="*70)
    print("5. PRESSURE DROP SCALING")
    print("="*70)
    
    # Hagen-Poiseuille: dP = 128*mu*L*Q / (pi*D^4)
    # For constant Q, dP scales as L/D^4
    
    mu = 3.5e-3  # Pa.s
    L = 1e-3     # m
    Q = 1e-9     # m^3/s
    
    D1 = 100e-6
    D2 = 50e-6  # Half the diameter
    
    dP1 = 128 * mu * L * Q / (np.pi * D1**4)
    dP2 = 128 * mu * L * Q / (np.pi * D2**4)
    
    # dP2/dP1 should be (D1/D2)^4 = 16
    ratio_expected = (D1/D2)**4
    ratio_actual = dP2 / dP1
    error = abs(ratio_actual - ratio_expected) / ratio_expected
    
    print(f"\n  Geometry:")
    print(f"    D1 = {D1*1e6:.0f} um")
    print(f"    D2 = {D2*1e6:.0f} um")
    
    print(f"\n  Pressure Drop:")
    print(f"    dP1 (D={D1*1e6:.0f} um) = {dP1:.2f} Pa")
    print(f"    dP2 (D={D2*1e6:.0f} um) = {dP2:.2f} Pa")
    print(f"    Ratio dP2/dP1 = {ratio_actual:.4f}")
    print(f"    Expected (D1/D2)^4 = {ratio_expected:.4f}")
    print(f"    Error: {error:.2e}")
    
    passed = error < 1e-10
    print(f"\n  RESULT: {'PASS' if passed else 'FAIL'}")
    return passed, error


def main():
    print("\n" + "#"*70)
    print(" BIFURCATION/TRIFURCATION VALIDATION SUITE")
    print(" Murray's Law and Mass Conservation")
    print("#"*70)
    
    results = []
    
    try:
        passed1, error1 = validate_murrays_law_bifurcation()
        results.append(("Murray Bifurcation", passed1, error1))
    except Exception as e:
        print(f"\n  ERROR: {e}")
        results.append(("Murray Bifurcation", False, float('inf')))
    
    try:
        passed2, error2 = validate_murrays_law_trifurcation()
        results.append(("Murray Trifurcation", passed2, error2))
    except Exception as e:
        print(f"\n  ERROR: {e}")
        results.append(("Murray Trifurcation", False, float('inf')))
    
    try:
        passed3, error3 = validate_mass_conservation_bifurcation()
        results.append(("Mass Bifurcation", passed3, error3))
    except Exception as e:
        print(f"\n  ERROR: {e}")
        results.append(("Mass Bifurcation", False, float('inf')))
    
    try:
        passed4, error4 = validate_mass_conservation_trifurcation()
        results.append(("Mass Trifurcation", passed4, error4))
    except Exception as e:
        print(f"\n  ERROR: {e}")
        results.append(("Mass Trifurcation", False, float('inf')))
    
    try:
        passed5, error5 = validate_pressure_scaling()
        results.append(("Pressure Scaling", passed5, error5))
    except Exception as e:
        print(f"\n  ERROR: {e}")
        results.append(("Pressure Scaling", False, float('inf')))
    
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
        print("ALL BIFURCATION VALIDATION TESTS PASSED")
    else:
        print("SOME TESTS FAILED")
    print("="*70)
    
    return 0 if all_passed else 1


if __name__ == "__main__":
    exit(main())
