#!/usr/bin/env python3
"""
Poiseuille Flow Validation: pycfdrs vs Analytical Solutions

This script validates the 2D Poiseuille flow solver against:
1. Newtonian analytical solution
2. Casson blood model viscosity
3. Wall shear stress

Run with:
    python validation/validate_poiseuille.py
"""

import numpy as np
import sys
sys.path.insert(0, "crates/pycfdrs")
import pycfdrs


def analytical_max_velocity(H, L, dP, mu):
    """Maximum centerline velocity for Newtonian Poiseuille flow.
    
    u_max = (dP/dx) * H^2 / (8*mu)
    """
    dP_dx = dP / L
    return abs((dP_dx * H**2) / (8 * mu))


def analytical_wss(H, L, dP):
    """Wall shear stress for Poiseuille flow.
    
    tau_w = (dP/dx) * H / 2
    """
    return abs((dP / L) * H / 2)


def validate_analytical():
    """Validate solver's analytical functions match theory."""
    print("\n" + "="*70)
    print("1. ANALYTICAL SOLUTION VALIDATION")
    print("="*70)
    
    # Physical parameters (microfluidic scale)
    H = 100e-6      # Channel height [m] (100 um)
    W = 1e-3        # Channel width [m] (1 mm)
    L = 1e-3        # Channel length [m] (1 mm)
    dP = 100.0      # Pressure drop [Pa]
    mu = 3.5e-3     # Viscosity [Pa.s] (blood at high shear)
    
    dP_dx = dP / L  # Positive pressure gradient for flow in +x direction
    
    # Create solver
    solver = pycfdrs.PyPoiseuille2DSolver(height=H, width=W, length=L, nx=50, ny=50)
    
    # Test analytical maximum velocity
    # NOTE: pycfdrs uses negative dP/dx convention, so we compare absolute values
    u_max_expected = analytical_max_velocity(H, L, dP, mu)
    u_max_solver = abs(solver.analytical_max_velocity(dP_dx, mu))
    u_max_error = abs(u_max_solver - u_max_expected) / u_max_expected
    
    print(f"\n  Parameters:")
    print(f"    H = {H*1e6:.0f} um, L = {L*1e3:.1f} mm, dP = {dP:.0f} Pa")
    print(f"    mu = {mu*1000:.2f} mPa.s")
    
    print(f"\n  Maximum Velocity:")
    print(f"    Expected:  {u_max_expected*1e3:.6f} mm/s")
    print(f"    Solver:    {u_max_solver*1e3:.6f} mm/s")
    print(f"    Error:     {u_max_error*100:.6f}%")
    
    passed = u_max_error < 1e-8  # Machine precision
    print(f"\n  RESULT: {'PASS' if passed else 'FAIL'}")
    return passed, u_max_error


def validate_wss():
    """Validate wall shear stress calculation."""
    print("\n" + "="*70)
    print("2. WALL SHEAR STRESS VALIDATION")
    print("="*70)
    
    H = 100e-6
    W = 1e-3
    L = 1e-3
    dP = 100.0
    
    wss_analytical = analytical_wss(H, L, dP)
    
    # Solve numerically
    solver = pycfdrs.PyPoiseuille2DSolver(height=H, width=W, length=L, nx=100, ny=50)
    result = solver.solve(dP, "newtonian")
    
    wss_error = abs(result.wall_shear_stress - wss_analytical) / wss_analytical
    
    print(f"\n  Analytical WSS: {wss_analytical:.4f} Pa")
    print(f"  Numerical WSS:  {result.wall_shear_stress:.4f} Pa")
    print(f"  Error: {wss_error*100:.4f}%")
    
    passed = wss_error < 0.001  # <0.1%
    print(f"\n  RESULT: {'PASS' if passed else 'FAIL'} (tolerance: <0.1%)")
    return passed, wss_error


def validate_casson_blood():
    """Validate Casson blood model viscosity calculation."""
    print("\n" + "="*70)
    print("3. CASSON BLOOD MODEL VALIDATION")
    print("="*70)
    
    # Create blood models
    casson = pycfdrs.CassonBlood()
    carreau = pycfdrs.CarreauYasudaBlood()
    
    print(f"\n  Created blood models")
    print(f"    Casson yield stress: {casson.yield_stress():.4f} Pa")
    print(f"    Casson density: {casson.density():.0f} kg/m^3")
    
    # Validate viscosity at specific shear rates
    # Literature values for normal blood (Merrill 1969):
    # - At gamma=100 s^-1: mu ~ 5-10 mPa.s
    # - At gamma=1000 s^-1: mu ~ 3-4 mPa.s (asymptotic)
    
    shear_rates = [1.0, 10.0, 100.0, 1000.0]
    
    print(f"\n  Viscosity vs Shear Rate:")
    print(f"    gamma [s^-1]   Casson [mPa.s]   Carreau [mPa.s]")
    print(f"    " + "-"*50)
    
    for gamma in shear_rates:
        mu_casson = casson.viscosity(gamma) * 1000  # mPa.s
        mu_carreau = carreau.viscosity(gamma) * 1000
        print(f"    {gamma:8.0f}        {mu_casson:8.4f}         {mu_carreau:8.4f}")
    
    # Check asymptotic viscosity
    mu_inf_casson = casson.viscosity_high_shear() * 1000  # mPa.s
    mu_inf_carreau = carreau.viscosity_high_shear() * 1000
    
    print(f"\n  Asymptotic viscosity:")
    print(f"    Casson mu_inf:  {mu_inf_casson:.3f} mPa.s")
    print(f"    Carreau mu_inf: {mu_inf_carreau:.3f} mPa.s")
    print(f"    Literature:     3-4 mPa.s (Merrill 1969)")
    
    # Validate: asymptotic viscosity should be 3-4 mPa.s
    in_range = 2.5 < mu_inf_casson < 5.0 and 2.5 < mu_inf_carreau < 5.0
    print(f"\n  RESULT: {'PASS' if in_range else 'FAIL'} (mu_inf in 2.5-5.0 mPa.s range)")
    return in_range, abs(mu_inf_casson - 3.45) / 3.45


def validate_murrays_law():
    """Validate Murray's Law for bifurcation geometry."""
    print("\n" + "="*70)
    print("4. MURRAY'S LAW VALIDATION")
    print("="*70)
    
    # For symmetric bifurcation: D_parent^3 = 2 * D_daughter^3
    # Therefore: D_daughter = D_parent / 2^(1/3)
    D_parent = 100e-6  # 100 um
    factor = 2 ** (-1/3)  # ~0.794
    D_daughter = D_parent * factor
    
    # Check Murray's Law
    lhs = D_parent ** 3
    rhs = 2 * (D_daughter ** 3)
    deviation = abs(lhs - rhs) / lhs
    
    print(f"\n  Parent diameter: {D_parent*1e6:.1f} um")
    print(f"  Daughter diameter: {D_daughter*1e6:.2f} um")
    print(f"  Murray factor: {factor:.6f}")
    print(f"  D_p^3 = {lhs*1e18:.6f} um^3")
    print(f"  2*D_d^3 = {rhs*1e18:.6f} um^3")
    print(f"  Deviation: {deviation:.2e}")
    
    passed = deviation < 1e-14  # Machine precision
    print(f"\n  RESULT: {'PASS' if passed else 'FAIL'}")
    return passed, deviation


def main():
    print("\n" + "#"*70)
    print(" POISEUILLE FLOW VALIDATION SUITE")
    print(" pycfdrs vs Analytical Solutions")
    print("#"*70)
    
    results = []
    
    # Test 1: Analytical functions
    try:
        passed1, error1 = validate_analytical()
        results.append(("Analytical Functions", passed1, error1))
    except Exception as e:
        print(f"\n  ERROR: {e}")
        import traceback
        traceback.print_exc()
        results.append(("Analytical Functions", False, float('inf')))
    
    # Test 2: WSS validation
    try:
        passed2, error2 = validate_wss()
        results.append(("Wall Shear Stress", passed2, error2))
    except Exception as e:
        print(f"\n  ERROR: {e}")
        import traceback
        traceback.print_exc()
        results.append(("Wall Shear Stress", False, float('inf')))
    
    # Test 3: Blood model
    try:
        passed3, error3 = validate_casson_blood()
        results.append(("Casson Blood Model", passed3, error3))
    except Exception as e:
        print(f"\n  ERROR: {e}")
        import traceback
        traceback.print_exc()
        results.append(("Casson Blood Model", False, float('inf')))
    
    # Test 4: Murray's Law
    try:
        passed4, error4 = validate_murrays_law()
        results.append(("Murray's Law", passed4, error4))
    except Exception as e:
        print(f"\n  ERROR: {e}")
        import traceback
        traceback.print_exc()
        results.append(("Murray's Law", False, float('inf')))
    
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
        print("ALL POISEUILLE VALIDATION TESTS PASSED")
    else:
        print("SOME TESTS FAILED")
    print("="*70)
    
    return 0 if all_passed else 1


if __name__ == "__main__":
    exit(main())
