#!/usr/bin/env python3
"""
Comprehensive cfd_python validation against analytical identities and invariants.

Validated components:
1. Blood rheology (Carreau-Yasuda and Casson)
2. Poiseuille 3D analytical consistency
3. Poiseuille 2D analytical consistency
4. Reynolds-number monotonicity under pressure scaling
5. Bifurcation mass conservation
6. Serpentine pressure scaling with segment count
"""

import math
import sys

try:
    import cfd_python
except ImportError:
    print("[FAIL] cfd_python not available. Run:")
    print("       maturin build --release --manifest-path crates/cfd-python/Cargo.toml")
    print("       python -m pip install --upgrade --force-reinstall target/wheels/cfd_python-*.whl")
    sys.exit(1)

MU_NEWTONIAN = 0.0035  # Pa*s (cfd_python newtonian convention)


def test_header(name: str) -> None:
    print("\n" + "=" * 80)
    print(f"TEST: {name}")
    print("=" * 80)


def pass_fail(passed: bool, message: str) -> bool:
    status = "[PASS]" if passed else "[FAIL]"
    print(f"{status} {message}")
    return passed


def rel_err(a: float, b: float) -> float:
    denom = max(abs(a), abs(b), 1.0e-30)
    return abs(a - b) / denom


results = {}

# ============================================================================
# TEST 1: Blood rheology invariants
# ============================================================================
test_header("Blood rheology invariants")

carreau = cfd_python.CarreauYasudaBlood()
casson = cfd_python.CassonBlood()

shear_rates = [10.0, 100.0, 500.0, 1000.0, 2000.0, 5000.0]
carreau_vals = [carreau.apparent_viscosity(g) for g in shear_rates]
casson_vals = [casson.apparent_viscosity(g) for g in shear_rates]

print(f"{'gamma_dot (s^-1)':>16} {'mu_carreau (mPa*s)':>22} {'mu_casson (mPa*s)':>20}")
print("-" * 62)
for g, mu_c, mu_s in zip(shear_rates, carreau_vals, casson_vals, strict=True):
    print(f"{g:16.1f} {mu_c * 1000.0:22.6f} {mu_s * 1000.0:20.6f}")

test1a = pass_fail(all(v > 0.0 for v in carreau_vals + casson_vals), "All apparent viscosities are positive")
test1b = pass_fail(carreau_vals[0] > carreau_vals[-1], "Carreau-Yasuda is shear-thinning")
test1c = pass_fail(casson.yield_stress() > 0.0, "Casson has positive yield stress")

results["1a. positive_viscosity"] = test1a
results["1b. carreau_shear_thinning"] = test1b
results["1c. casson_yield_stress"] = test1c

# ============================================================================
# TEST 2: Poiseuille 3D analytical consistency
# ============================================================================
test_header("Poiseuille 3D analytical consistency")

D = 100e-6
L = 10e-3
pressure_drop = -4000.0  # P_out - P_in; negative => forward flow (Q > 0)

solver_3d = cfd_python.Poiseuille3DSolver(D, L, 20, 16, 30)
res_3d = solver_3d.solve(pressure_drop, "newtonian")

q_analytic_3d = (math.pi * D**4 * abs(pressure_drop)) / (128.0 * MU_NEWTONIAN * L)
q_error_3d = rel_err(abs(res_3d.flow_rate), q_analytic_3d)

print(f"flow_rate_numeric = {res_3d.flow_rate:.12e} m^3/s")
print(f"flow_rate_analytic = {q_analytic_3d:.12e} m^3/s")
print(f"relative_error     = {q_error_3d:.3e}")

test2a = pass_fail(q_error_3d < 1.0e-12, "3D flow-rate matches Hagen-Poiseuille")
test2b = pass_fail(res_3d.flow_rate > 0.0, "3D sign convention: pressure_drop < 0 gives forward flow")

results["2a. poiseuille_3d_flow_rate"] = test2a
results["2b. poiseuille_3d_sign"] = test2b

# ============================================================================
# TEST 3: Poiseuille 2D analytical consistency
# ============================================================================
test_header("Poiseuille 2D analytical consistency")

height = 100e-6
width = 100e-6
solver_2d = cfd_python.Poiseuille2DSolver(height, width, L, 60, 60)
res_2d = solver_2d.solve(pressure_drop, "newtonian")

pressure_gradient = pressure_drop / L
q_analytic_2d = solver_2d.analytical_flow_rate(pressure_gradient, MU_NEWTONIAN)
q_error_2d = rel_err(res_2d.flow_rate, q_analytic_2d)

print(f"flow_rate_numeric = {res_2d.flow_rate:.12e} m^3/s")
print(f"flow_rate_analytic = {q_analytic_2d:.12e} m^3/s")
print(f"relative_error     = {q_error_2d:.3e}")
print(f"Reynolds_number    = {res_2d.reynolds_number:.6f}")

test3a = pass_fail(q_error_2d < 1.0e-12, "2D flow-rate matches analytical parallel-plate solution")
test3b = pass_fail(res_2d.reynolds_number > 0.0, "2D Reynolds number is physically positive")

results["3a. poiseuille_2d_flow_rate"] = test3a
results["3b. poiseuille_2d_reynolds"] = test3b

# ============================================================================
# TEST 4: Reynolds scaling monotonicity
# ============================================================================
test_header("Reynolds scaling monotonicity")

dp_values = [-1000.0, -2000.0, -4000.0, -8000.0]
re_values = []
for dp in dp_values:
    r = solver_3d.solve(dp, "newtonian")
    re_values.append(abs(r.reynolds_number))
    print(f"dp={dp:8.1f} Pa -> |Re|={abs(r.reynolds_number):.8f}")

monotonic = all(re_values[i + 1] > re_values[i] for i in range(len(re_values) - 1))
ratio_1 = re_values[1] / re_values[0]
ratio_2 = re_values[2] / re_values[1]
ratio_3 = re_values[3] / re_values[2]
linear_like = all(abs(r - 2.0) < 0.02 for r in (ratio_1, ratio_2, ratio_3))

print(f"ratio(Re2/Re1)={ratio_1:.6f}, ratio(Re3/Re2)={ratio_2:.6f}, ratio(Re4/Re3)={ratio_3:.6f}")

test4a = pass_fail(monotonic, "Reynolds magnitude increases monotonically with |pressure_drop|")
test4b = pass_fail(linear_like, "Reynolds scaling is approximately linear under pressure doubling")

results["4a. reynolds_monotonic"] = test4a
results["4b. reynolds_linear_scaling"] = test4b

# ============================================================================
# TEST 5: Bifurcation mass conservation
# ============================================================================
test_header("Bifurcation mass conservation")

bif = cfd_python.BifurcationSolver(200e-6, 140e-6, 140e-6, 5e-3)
q_in = 2e-9
bif_res = bif.solve(q_in, 40.0, "casson")

q_out = bif_res.q_1 + bif_res.q_2
mass_error = rel_err(q_out, bif_res.q_parent)
split = bif_res.q_1 / max(q_out, 1.0e-30)

print(f"q_parent = {bif_res.q_parent:.12e} m^3/s")
print(f"q_1      = {bif_res.q_1:.12e} m^3/s")
print(f"q_2      = {bif_res.q_2:.12e} m^3/s")
print(f"mass_err = {mass_error:.3e}")
print(f"split_1  = {split:.6f}")

test5a = pass_fail(mass_error < 1.0e-12, "Bifurcation conserves mass")
test5b = pass_fail(abs(split - 0.5) < 1.0e-12, "Symmetric bifurcation yields equal daughter split")

results["5a. bifurcation_mass_conservation"] = test5a
results["5b. bifurcation_equal_split"] = test5b

# ============================================================================
# TEST 6: Serpentine pressure scaling
# ============================================================================
test_header("Serpentine pressure scaling")

velocity_inlet = 0.05
segment_counts = [5, 10, 20]
pressure_drops = []

for nseg in segment_counts:
    serp = cfd_python.SerpentineSolver1D(100e-6, 100e-6, 500e-6, nseg, 200e-6)
    serp_res = serp.solve(velocity_inlet, "casson")
    pressure_drops.append(serp_res.pressure_drop)
    print(f"segments={nseg:2d} -> dP={serp_res.pressure_drop:.6f} Pa")

monotonic_dp = all(pressure_drops[i + 1] > pressure_drops[i] for i in range(len(pressure_drops) - 1))
dp_per_segment = [dp / nseg for dp, nseg in zip(pressure_drops, segment_counts, strict=True)]
spread = (max(dp_per_segment) - min(dp_per_segment)) / max(dp_per_segment)

print(
    "dP/segment = "
    + ", ".join(f"{val:.6f}" for val in dp_per_segment)
    + f" Pa (spread={spread*100:.2f}%)"
)

test6a = pass_fail(monotonic_dp, "Serpentine pressure drop increases with segment count")
test6b = pass_fail(spread < 0.05, "Pressure drop per segment remains near-constant (<5% spread)")

results["6a. serpentine_monotonic_dp"] = test6a
results["6b. serpentine_linear_dp_per_segment"] = test6b

# ============================================================================
# SUMMARY
# ============================================================================
print("\n" + "=" * 80)
print("VALIDATION SUMMARY")
print("=" * 80)

passed = sum(1 for v in results.values() if v)
total = len(results)

for name, ok in results.items():
    print(f"{'[PASS]' if ok else '[FAIL]'} {name}")

print(f"\n{passed}/{total} checks passed ({(100.0 * passed / total):.1f}%)")

if passed == total:
    print("[PASS] Comprehensive cfd_python validation succeeded.")
    sys.exit(0)

print("[FAIL] Comprehensive cfd_python validation has failing checks.")
sys.exit(1)
