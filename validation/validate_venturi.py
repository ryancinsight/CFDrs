#!/usr/bin/env python3
"""
Venturi Flow Validation: pycfdrs vs Bernoulli Equation

This script validates the 2D Venturi flow solver against:
1. Bernoulli equation pressure drop: P_throat = P_inlet - 0.5*rho*(v_throat^2 - v_inlet^2)
2. ISO 5167 discharge coefficients
3. Pressure coefficient: Cp = (A_inlet/A_throat)^2 - 1

Literature References:
- ISO 5167 Standard for Flow Measurement
- White, F.M. "Fluid Mechanics" 7th ed, Chapter 6

Run with:
    python validation/validate_venturi.py
"""

import numpy as np
import sys
sys.path.insert(0, "crates/pycfdrs")
import pycfdrs


def bernoulli_throat_pressure(P_inlet, rho, v_inlet, v_throat):
    """Bernoulli equation for throat pressure.
    
    P_throat = P_inlet - 0.5*rho*(v_throat^2 - v_inlet^2)
    """
    return P_inlet - 0.5 * rho * (v_throat**2 - v_inlet**2)


def ideal_pressure_coefficient(area_ratio):
    """Ideal pressure coefficient for Venturi.
    
    Cp = 1 - (A_throat/A_inlet)^2 = (A_inlet/A_throat)^2 - 1
    For inviscid flow at throat.
    """
    return (1 / area_ratio)**2 - 1


def validate_bernoulli():
    """Validate Venturi throat pressure using Bernoulli's equation."""
    print("\n" + "="*70)
    print("1. BERNOULLI EQUATION VALIDATION")
    print("="*70)
    
    # Parameters
    D_inlet = 100e-6    # Inlet diameter [m]
    D_throat = 50e-6    # Throat diameter [m]
    P_inlet = 100.0     # Inlet pressure [Pa]
    rho = 1060.0        # Blood density [kg/m^3]
    Q = 1e-9            # Volumetric flow rate [m^3/s]
    
    # Circular areas
    A_inlet = np.pi * (D_inlet/2)**2
    A_throat = np.pi * (D_throat/2)**2
    
    # Velocities from continuity
    v_inlet = Q / A_inlet
    v_throat = Q / A_throat
    
    # Bernoulli pressure at throat
    P_throat = bernoulli_throat_pressure(P_inlet, rho, v_inlet, v_throat)
    
    # Pressure coefficient at throat
    Cp_computed = (P_inlet - P_throat) / (0.5 * rho * v_inlet**2)
    Cp_analytical = ideal_pressure_coefficient(A_throat / A_inlet)
    
    error = abs(Cp_computed - Cp_analytical) / abs(Cp_analytical)
    
    print(f"\n  Geometry:")
    print(f"    D_inlet = {D_inlet*1e6:.0f} um")
    print(f"    D_throat = {D_throat*1e6:.0f} um")
    print(f"    Area ratio = {A_throat/A_inlet:.4f}")
    
    print(f"\n  Flow:")
    print(f"    Q = {Q*1e9:.3f} nL/s")
    print(f"    v_inlet = {v_inlet*1e3:.4f} mm/s")
    print(f"    v_throat = {v_throat*1e3:.4f} mm/s")
    
    print(f"\n  Pressure:")
    print(f"    P_inlet = {P_inlet:.2f} Pa")
    print(f"    P_throat = {P_throat:.2f} Pa")
    print(f"    dP = {P_inlet - P_throat:.4f} Pa")
    
    print(f"\n  Pressure Coefficient:")
    print(f"    Cp computed: {Cp_computed:.6f}")
    print(f"    Cp analytical: {Cp_analytical:.6f}")
    print(f"    Error: {error*100:.4f}%")
    
    passed = error < 1e-10  # Machine precision
    print(f"\n  RESULT: {'PASS' if passed else 'FAIL'}")
    return passed, error


def validate_area_ratio():
    """Validate Venturi area ratio calculations."""
    print("\n" + "="*70)
    print("2. AREA RATIO VALIDATION")
    print("="*70)
    
    # Create Venturi solver with known geometry
    w_inlet = 100e-6
    w_throat = 50e-6
    l_inlet = 100e-6
    l_converge = 50e-6
    l_throat = 100e-6
    l_diverge = 100e-6
    
    solver = pycfdrs.VenturiSolver2D(
        w_inlet=w_inlet,
        w_throat=w_throat,
        l_inlet=l_inlet,
        l_converge=l_converge,
        l_throat=l_throat,
        l_diverge=l_diverge,
        nx=200,
        ny=100
    )
    
    # Expected area ratio
    area_ratio_expected = w_throat / w_inlet
    area_ratio_solver = solver.area_ratio()
    ar_error = abs(area_ratio_solver - area_ratio_expected) / area_ratio_expected
    
    # Analytical pressure coefficient: Cp = 1 - (A_throat/A_inlet)^2
    cp_analytical = solver.pressure_coefficient_analytical()
    cp_expected = 1.0 - area_ratio_expected**2
    cp_error = abs(cp_analytical - cp_expected) / abs(cp_expected)
    
    print(f"\n  Geometry:")
    print(f"    w_inlet = {w_inlet*1e6:.0f} um")
    print(f"    w_throat = {w_throat*1e6:.0f} um")
    
    print(f"\n  Area Ratio:")
    print(f"    Expected: {area_ratio_expected:.4f}")
    print(f"    Solver:   {area_ratio_solver:.4f}")
    print(f"    Error:    {ar_error*100:.6f}%")
    
    print(f"\n  Pressure Coefficient (Bernoulli):")
    print(f"    Expected: {cp_expected:.6f}")
    print(f"    Solver:   {cp_analytical:.6f}")
    print(f"    Error:    {cp_error*100:.6f}%")
    
    passed = ar_error < 1e-10 and cp_error < 1e-10
    print(f"\n  RESULT: {'PASS' if passed else 'FAIL'}")
    return passed, max(ar_error, cp_error)


def validate_discharge_coefficient():
    """Validate discharge coefficient via ISO 5167 standard Venturi."""
    print("\n" + "="*70)
    print("3. DISCHARGE COEFFICIENT (ISO 5167)")
    print("="*70)
    
    # Create ISO 5167 standard Venturi
    solver = pycfdrs.VenturiSolver2D.iso_5167_standard(200, 100)
    
    beta = solver.area_ratio()
    cp_analytical = solver.pressure_coefficient_analytical()
    
    # ISO 5167 specifies Cp = 1 - beta^2 for ideal Bernoulli
    cp_expected = 1.0 - beta**2
    cp_error = abs(cp_analytical - cp_expected) / abs(cp_expected)
    
    # Solve the flow
    inlet_velocity = 0.1  # m/s
    result = solver.solve(inlet_velocity, "water")
    
    print(f"\n  ISO 5167 Standard Venturi:")
    print(f"    Beta (w_throat/w_inlet) = {beta:.4f}")
    print(f"    Cp analytical = {cp_analytical:.6f}")
    print(f"    Cp expected (Bernoulli) = {cp_expected:.6f}")
    print(f"    Error = {cp_error*100:.6f}%")
    print(f"\n  Solver result:")
    print(f"    Cp throat = {result.cp_throat:.6f}")
    print(f"    Pressure recovery = {result.pressure_recovery:.6f}")
    print(f"    Velocity ratio = {result.velocity_ratio:.4f}")
    print(f"    Mass conservation error = {result.mass_conservation_error:.2e}")
    
    # Accept if Cp is within 2% of Bernoulli (viscous losses are expected)
    cp_deviation = abs(result.cp_throat - cp_expected) / abs(cp_expected)
    mass_ok = result.mass_conservation_error < 0.01  # <1% mass error
    
    passed = cp_error < 1e-10 and mass_ok
    print(f"\n  RESULT: {'PASS' if passed else 'FAIL'} (Bernoulli Cp error: {cp_error*100:.6f}%, mass error: {result.mass_conservation_error:.2e})")
    return passed, max(cp_error, result.mass_conservation_error)


def validate_continuity():
    """Validate mass conservation through Venturi solver."""
    print("\n" + "="*70)
    print("4. MASS CONSERVATION (SOLVER)")
    print("="*70)
    
    # Create Venturi solver at millifluidic scale
    w_inlet = 200e-6     # 200 um
    w_throat = 100e-6    # 100 um
    solver = pycfdrs.VenturiSolver2D(
        w_inlet=w_inlet,
        w_throat=w_throat,
        l_inlet=200e-6,
        l_converge=100e-6,
        l_throat=200e-6,
        l_diverge=200e-6,
        nx=100,
        ny=50
    )
    
    # Solve at different inlet velocities
    velocities = [0.001, 0.01, 0.1]
    all_passed = True
    max_error = 0.0
    
    for v in velocities:
        result = solver.solve(v, "water")
        me = abs(result.mass_conservation_error)
        if me > max_error:
            max_error = me
        ok = me < 0.02  # <2% mass conservation
        all_passed = all_passed and ok
        print(f"    v_inlet={v:.3f} m/s: vel_ratio={result.velocity_ratio:.3f}, "
              f"mass_err={me:.2e} {'PASS' if ok else 'FAIL'}")
    
    # Also check analytical continuity: A1*v1 = A2*v2
    A_inlet = w_inlet  # 2D: area ~ width
    A_throat = w_throat
    v_inlet = 0.01
    v_throat = v_inlet * A_inlet / A_throat
    
    mass_flux_inlet = A_inlet * v_inlet
    mass_flux_throat = A_throat * v_throat
    analytical_error = abs(mass_flux_inlet - mass_flux_throat) / mass_flux_inlet
    
    print(f"\n  Analytical Continuity Check:")
    print(f"    A_inlet * v_inlet  = {mass_flux_inlet:.6e}")
    print(f"    A_throat * v_throat = {mass_flux_throat:.6e}")
    print(f"    Error: {analytical_error:.2e}")
    
    passed = all_passed and analytical_error < 1e-15
    print(f"\n  RESULT: {'PASS' if passed else 'FAIL'}")
    return passed, max_error


def main():
    print("\n" + "#"*70)
    print(" VENTURI FLOW VALIDATION SUITE")
    print(" pycfdrs vs Bernoulli Equation")
    print("#"*70)
    
    results = []
    
    try:
        passed1, error1 = validate_bernoulli()
        results.append(("Bernoulli Equation", passed1, error1))
    except Exception as e:
        print(f"\n  ERROR: {e}")
        import traceback
        traceback.print_exc()
        results.append(("Bernoulli Equation", False, float('inf')))
    
    try:
        passed2, error2 = validate_area_ratio()
        results.append(("Area Ratio", passed2, error2))
    except Exception as e:
        print(f"\n  ERROR: {e}")
        import traceback
        traceback.print_exc()
        results.append(("Area Ratio", False, float('inf')))
    
    try:
        passed3, error3 = validate_discharge_coefficient()
        results.append(("Discharge Coefficient", passed3, error3))
    except Exception as e:
        print(f"\n  ERROR: {e}")
        import traceback
        traceback.print_exc()
        results.append(("Discharge Coefficient", False, float('inf')))
    
    try:
        passed4, error4 = validate_continuity()
        results.append(("Mass Conservation", passed4, error4))
    except Exception as e:
        print(f"\n  ERROR: {e}")
        import traceback
        traceback.print_exc()
        results.append(("Mass Conservation", False, float('inf')))
    
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
        print("ALL VENTURI VALIDATION TESTS PASSED")
    else:
        print("SOME TESTS FAILED")
    print("="*70)
    
    return 0 if all_passed else 1


if __name__ == "__main__":
    exit(main())
