# CAVITY SIMPLEC SOLVER - DEBUG STATE AND NEXT STEPS
# Last updated: Current session

## WORKSPACE
- Path: c:\Users\RyanClanton\gcli\software\millifluidic_design\CFDrs
- Python: .\.venv\Scripts\python.exe
- Build: maturin develop --manifest-path crates/pycfdrs/Cargo.toml (add --release for opt builds)

## ALL OTHER VALIDATIONS PASS
- external_reference_comparison.py: 6/6 PASSED
- comprehensive_validation.py: 31/31 PASSED  
- All individual scripts: PASSED
- Rust validation: 50/50 PASSED

## COMPLETED FIXES THIS SESSION
1. Rhie-Chow coefficient fix in solver.rs (3 places) - changed ap_consistent -> full ap for Rhie-Chow
2. Source term overwrite fix in coefficients.rs - changed *= to += for deferred correction
3. Pre-existing venturi_flow.rs compilation fix

## CAVITY STATUS: STILL NOT CONVERGING
- 9x9 grid: residual goes from 0.64 to 3.42 and plateaus over 2000 steps
- 17x17 grid: residual stuck at ~4.53 from step 0

## ROOT CAUSE ANALYSIS - ADDITIONAL ISSUES FOUND

### Issue 1: Rhie-Chow Transient Correction (rhie_chow.rs)
The transient correction formula is:
```rust
let transient_factor = dt / two;
u_f += transient_factor * (u_bar - u_bar_old);
```
This should probably be: d_face * rho/dt * (u_bar_old - u_bar) or just removed for simpler debugging.
The sign (u_bar - u_bar_old vs u_bar_old - u_bar) and the factor (dt/2 vs d_face*rho/dt) may be wrong.

### Issue 2: Pressure Equation vs Velocity Correction Mismatch
- Pressure equation uses FACE-BASED d coefficients from Rhie-Chow (harmonic avg of V/ap_full)
- Velocity correction uses V/ap_consistent with CENTERED difference for dp'/dx
- These two formulations should be consistent for mathematical convergence

### Issue 3: Debug Build Performance
- Need --release flag for maturin build; debug is 10-100x slower
- Command: maturin develop --release --manifest-path crates/pycfdrs/Cargo.toml

## RECOMMENDED APPROACH: Python Reference Solver
Instead of continuing to debug the complex Rust SIMPLEC in blind, write a standalone
Python SIMPLE solver (50-100 lines with NumPy) for the lid-driven cavity problem.
This provides:
1. A known-working reference to compare against
2. Ability to print intermediate values for debugging
3. Step-by-step comparison of what pycfdrs computes vs what Python computes

### Python Reference Solver Plan:
- Grid: 33x33 (matching pycfdrs)
- Algorithm: Basic SIMPLE (not SIMPLEC) with staggered grid layout
- BCs: lid velocity 1.0 at top, no-slip elsewhere
- Re=100 (rho=1, mu=0.01, L=1, U=1)
- Compare centerline velocity profiles with Ghia et al. (1982)
- If Python SIMPLE converges correctly, compare step-by-step with pycfdrs

### Alternative: Disable Rhie-Chow and Use Staggered-like Approach
The Rhie-Chow interpolation adds complexity. Could try:
1. Set config.use_rhie_chow = false in PyCavitySolver2D
2. Use the simpler pressure correction without face velocities
3. This bypasses all Rhie-Chow issues and tests the basic SIMPLE algorithm

## KEY SOURCE FILES
- crates/cfd-2d/src/simplec_pimple/solver.rs (1011 lines) - SIMPLEC/PIMPLE solver
- crates/cfd-2d/src/physics/momentum/coefficients.rs - FVM momentum discretization
- crates/cfd-2d/src/physics/momentum/solver.rs (570 lines) - MomentumSolver
- crates/cfd-2d/src/pressure_velocity/pressure.rs (608 lines) - Pressure correction
- crates/cfd-2d/src/pressure_velocity/rhie_chow.rs - Rhie-Chow interpolation
- crates/pycfdrs/src/solver_2d.rs (740 lines) - Python bindings

## GHIA BENCHMARK DATA (Re=100)
# U-velocity along vertical centerline at x=0.5:
# y=0.0000: u=0.0000 (no-slip)
# y=0.0547: u=-0.0372
# y=0.0625: u=-0.0419
# y=0.0703: u=-0.0477
# y=0.1016: u=-0.0643
# y=0.1719: u=-0.1015
# y=0.2813: u=-0.1566
# y=0.4531: u=-0.2109
# y=0.5000: u=-0.2058
# y=0.6172: u=-0.1364
# y=0.7344: u=0.0033
# y=0.8516: u=0.2315
# y=0.9531: u=0.6872
# y=0.9609: u=0.7372
# y=0.9688: u=0.7887
# y=1.0000: u=1.0000 (lid)
