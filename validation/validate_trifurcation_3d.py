from pycfdrs import Trifurcation3DSolver

def validate_trifurcation():
    print("\n" + "="*70)
    print("3D Trifurcation Validation: Newtonian Balance")
    print("="*70 + "\n")
    
    # 2. Config
    # Equal split geometry: Q_in = Q_d1 + Q_d2 + Q_d3
    q_total = 1e-6
    d_parent = 1e-3
    d_daughter = 0.5e-3
    
    # 3. Solve
    # Increase resolution slightly for reliability
    solver = Trifurcation3DSolver(d_parent, d_daughter, length=5e-3)
    result = solver.solve(flow_rate=q_total, blood_type="newtonian")
    
    # result = solver.solve(...)
    print(f"Flow Rates [m3/s]:")
    print(f"  Inlet:      {result.flow_rates[0]:.2e}")
    print(f"  Daughter 1: {result.flow_rates[1]:.2e}")
    print(f"  Daughter 2: {result.flow_rates[2]:.2e}")
    print(f"  Daughter 3: {result.flow_rates[3]:.2e}")
    
    # 4. Check Mass Balance
    q_sum_out = sum(result.flow_rates[1:])
    mass_error = abs(result.flow_rates[0] - q_sum_out) / result.flow_rates[0]
    
    print(f"\nMass Balance Error: {mass_error:.2%}")
    
    # 5. Check Symmetry (Should be exactly 1/3 if symmetric)
    q_ideal = q_total / 3.0
    symmetry_error = max(abs(q - q_ideal) for q in result.flow_rates[1:]) / q_ideal
    print(f"Symmetry Error:     {symmetry_error:.2%}")

    threshold = 0.05
    passed = (mass_error < 1e-6) and (symmetry_error < threshold)
    
    print(f"\nValidation Status: {'PASSED' if passed else 'FAILED'}")
    
    return passed

if __name__ == "__main__":
    validate_trifurcation()
