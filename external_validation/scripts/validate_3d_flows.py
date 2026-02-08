#!/usr/bin/env python3
"""
3D Flow Validation

Validates CFD-rs 3D solvers against:
1. Analytical solutions (Hagen-Poiseuille)
2. Murray's law for bifurcations
3. Wall shear stress correlations

Blood flow specific validations:
- Casson and Carreau-Yasuda rheology
- Fåhræus-Lindqvist effect in microvessels
- Physiological wall shear stress ranges
"""

import argparse
import json
import sys
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Dict, List
import numpy as np

try:
    import pycfdrs
    HAS_PYCFDRS = True
except ImportError:
    HAS_PYCFDRS = False
    print("Error: pycfdrs not available. Build with: maturin develop")
    sys.exit(1)


try:
    import matplotlib.pyplot as plt
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False


@dataclass
class ValidationResult3D:
    """Result of 3D flow validation"""
    test_name: str
    cfrs_value: float
    reference_value: float
    absolute_error: float
    relative_error: float
    tolerance: float
    passed: bool
    units: str
    description: str
    
    def to_dict(self) -> dict:
        return asdict(self)


class Validator3D:
    """
    Validates 3D CFD-rs solvers against analytical and literature solutions.
    """
    
    def __init__(self, output_dir: str = "external_validation/results"):
        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        self.results: List[ValidationResult3D] = []
        
    def validate_poiseuille_3d(self) -> ValidationResult3D:
        """
        Validate 3D Poiseuille flow against Hagen-Poiseuille analytical.
        
        Analytical solution for laminar pipe flow:
        u(r) = (ΔP/4μL)(R² - r²)
        Q = πR⁴ΔP / (8μL)
        """
        print("\nValidating 3D Poiseuille flow...")
        
        # Test parameters
        diameter = 100e-6  # 100 μm
        radius = diameter / 2
        length = 1e-3  # 1 mm
        viscosity = 0.001  # Pa·s (water)
        delta_p = 100.0  # Pa
        
        # Analytical solution
        q_analytical = np.pi * radius**4 * delta_p / (8 * viscosity * length)
        
        # CFD-rs solution
        solver = pycfdrs.Poiseuille3DSolver(
            diameter=diameter,
            length=length,
            nr=10,
            ntheta=16,
            nz=20
        )
        
        result = solver.solve(pressure_drop=delta_p, blood_type="casson")
        q_cfrs = result.flow_rate
        
        # Calculate error
        abs_error = abs(q_cfrs - q_analytical)
        rel_error = abs_error / q_analytical
        
        result_val = ValidationResult3D(
            test_name="Poiseuille 3D Flow Rate",
            cfrs_value=q_cfrs,
            reference_value=q_analytical,
            absolute_error=abs_error,
            relative_error=rel_error,
            tolerance=0.05,  # 5% tolerance for 3D
            passed=rel_error < 0.05,
            units="m³/s",
            description="3D pipe flow vs Hagen-Poiseuille analytical"
        )
        self.results.append(result_val)
        
        print(f"  Analytical Q: {q_analytical:.4e} m³/s")
        print(f"  CFD-rs Q:     {q_cfrs:.4e} m³/s")
        print(f"  Rel. error:   {rel_error*100:.2f}%")
        print(f"  Status:       {'PASS' if result_val.passed else 'FAIL'}")
        
        return result_val
    
    def validate_bifurcation_murray(self) -> ValidationResult3D:
        """
        Validate 3D bifurcation against Murray's law.
        
        Murray's law: D_parent³ = D_daughter1³ + D_daughter2³
        Optimal for minimizing power consumption.
        """
        print("\nValidating 3D bifurcation (Murray's law)...")
        
        # Murray's law: D_daughter = D_parent / 2^(1/3)
        d_parent = 100e-6
        d_daughter = d_parent / (2**(1/3))  # ~79.37 μm
        
        # Check Murray compliance
        lhs = d_parent**3
        rhs = 2 * d_daughter**3
        murray_deviation = abs(lhs - rhs) / lhs
        
        # CFD-rs simulation
        solver = pycfdrs.Bifurcation3DSolver(
            d_parent=d_parent,
            d_daughter1=d_daughter,
            d_daughter2=d_daughter,
            angle=45.0,
            length=1e-3,
            nx=30, ny=30, nz=30
        )
        
        result = solver.solve(flow_rate=1e-9, blood_type="casson")
        
        # Check mass conservation
        mass_error = result.mass_conservation_error
        
        # Combined validation
        passed = murray_deviation < 0.01 and mass_error < 1e-10
        
        result_val = ValidationResult3D(
            test_name="Bifurcation Murray's Law",
            cfrs_value=mass_error,
            reference_value=0.0,
            absolute_error=mass_error,
            relative_error=mass_error,
            tolerance=1e-10,
            passed=passed,
            units="dimensionless",
            description="Murray's law compliance and mass conservation"
        )
        self.results.append(result_val)
        
        print(f"  D_parent: {d_parent*1e6:.2f} μm")
        print(f"  D_daughter: {d_daughter*1e6:.2f} μm")
        print(f"  Murray deviation: {murray_deviation*100:.4f}%")
        print(f"  Mass error: {mass_error:.2e}")
        print(f"  Flow split: {result.flow_split_ratio:.4f}")
        print(f"  Status: {'PASS' if result_val.passed else 'FAIL'}")
        
        return result_val
    
    def validate_wall_shear_stress(self) -> ValidationResult3D:
        """
        Validate wall shear stress in 3D pipe flow.
        
        Analytical WSS for Poiseuille flow:
        τ_w = R * ΔP / (2L) = 4μQ / (πR³)
        
        Physiological range:
        - Capillaries: 1-5 Pa
        - Arteries: 0.5-3 Pa
        - Veins: 0.1-0.6 Pa
        """
        print("\nValidating wall shear stress...")
        
        # Parameters for microvessel
        diameter = 50e-6  # 50 μm (capillary)
        radius = diameter / 2
        length = 1e-3
        viscosity = 0.003  # Blood viscosity
        delta_p = 200.0  # Pa
        
        # Analytical WSS
        tau_w_analytical = radius * delta_p / (2 * length)
        
        # CFD-rs simulation
        solver = pycfdrs.Poiseuille3DSolver(
            diameter=diameter,
            length=length,
            nr=8,
            ntheta=12,
            nz=15
        )
        
        result = solver.solve(pressure_drop=delta_p, blood_type="casson")
        tau_w_cfrs = result.wall_shear_stress
        
        # Calculate error
        rel_error = abs(tau_w_cfrs - tau_w_analytical) / tau_w_analytical
        
        # Check physiological range (capillaries: 1-5 Pa)
        physiological = 1.0 <= tau_w_cfrs <= 5.0
        
        result_val = ValidationResult3D(
            test_name="Wall Shear Stress (3D)",
            cfrs_value=tau_w_cfrs,
            reference_value=tau_w_analytical,
            absolute_error=abs(tau_w_cfrs - tau_w_analytical),
            relative_error=rel_error,
            tolerance=0.10,  # 10% tolerance
            passed=rel_error < 0.10 and physiological,
            units="Pa",
            description="WSS in microvessel vs analytical"
        )
        self.results.append(result_val)
        
        print(f"  Analytical WSS: {tau_w_analytical:.4f} Pa")
        print(f"  CFD-rs WSS:     {tau_w_cfrs:.4f} Pa")
        print(f"  Rel. error:     {rel_error*100:.2f}%")
        print(f"  Physiological:  {physiological} (1-5 Pa for capillaries)")
        print(f"  Status:         {'PASS' if result_val.passed else 'FAIL'}")
        
        return result_val
    
    def validate_trifurcation_flow(self) -> ValidationResult3D:
        """
        Validate 3D trifurcation flow distribution.
        
        For symmetric trifurcation:
        Q_daughter = Q_parent / 3
        """
        print("\nValidating 3D trifurcation flow...")
        
        d_parent = 100e-6
        d_daughter = d_parent / (3**(1/3))  # ~69.3 μm
        q_parent = 1e-9  # 1 nL/s
        
        # Expected flow per daughter
        q_expected = q_parent / 3
        
        # CFD-rs simulation
        solver = pycfdrs.Trifurcation3DSolver(
            d_parent=d_parent,
            d_daughter=d_daughter,
            length=1e-3
        )
        
        result = solver.solve(flow_rate=q_parent, blood_type="casson")
        
        # Get daughter flows (indices 1, 2, 3)
        q_daughters = result.flow_rates[1:]
        q_mean = np.mean(q_daughters)
        
        # Calculate error
        rel_error = abs(q_mean - q_expected) / q_expected
        
        # Check conservation
        q_out = sum(q_daughters)
        mass_error = abs(q_parent - q_out) / q_parent
        
        result_val = ValidationResult3D(
            test_name="Trifurcation Flow Split",
            cfrs_value=q_mean,
            reference_value=q_expected,
            absolute_error=abs(q_mean - q_expected),
            relative_error=rel_error,
            tolerance=0.05,
            passed=rel_error < 0.05 and mass_error < 1e-10,
            units="m³/s",
            description="Symmetric trifurcation flow distribution"
        )
        self.results.append(result_val)
        
        print(f"  Q_parent: {q_parent:.4e} m³/s")
        print(f"  Q_expected per daughter: {q_expected:.4e} m³/s")
        print(f"  Q_daughters: {[f'{q:.4e}' for q in q_daughters]}")
        print(f"  Mean error: {rel_error*100:.2f}%")
        print(f"  Mass error: {mass_error:.2e}")
        print(f"  Status: {'PASS' if result_val.passed else 'FAIL'}")
        
        return result_val
    
    def validate_blood_rheology(self) -> ValidationResult3D:
        """
        Validate non-Newtonian blood rheology models.
        
        Casson model:
        μ_app = (√τ_y/√γ̇ + √μ_∞)²
        
        For normal blood:
        - τ_y = 0.0056 Pa (yield stress)
        - μ_∞ = 0.00345 Pa·s (infinite-shear viscosity)
        """
        print("\nValidating blood rheology models...")
        
        # Test shear rates
        shear_rates = [1.0, 10.0, 100.0, 1000.0]  # s⁻¹
        
        # Casson parameters
        tau_y = 0.0056  # Pa
        mu_inf = 0.00345  # Pa·s
        
        errors = []
        for gamma in shear_rates:
            # Analytical Casson viscosity
            sqrt_term = np.sqrt(tau_y / gamma) + np.sqrt(mu_inf)
            mu_analytical = sqrt_term**2
            
            # CFD-rs viscosity (through Poiseuille solver)
            solver = pycfdrs.Poiseuille3DSolver(
                diameter=100e-6,
                length=1e-3,
                nr=5, ntheta=8, nz=10
            )
            
            # Use Casson blood model
            result = solver.solve(pressure_drop=100.0, blood_type="casson")
            
            # Approximate viscosity from Reynolds number
            # Re = ρuD/μ → μ = ρuD/Re
            rho = 1060.0  # Blood density
            u = result.max_velocity
            D = 100e-6
            Re = result.reynolds_number
            mu_numerical = rho * u * D / Re
            
            errors.append(abs(mu_numerical - mu_analytical) / mu_analytical)
        
        max_error = max(errors)
        
        result_val = ValidationResult3D(
            test_name="Blood Rheology (Casson)",
            cfrs_value=max_error,
            reference_value=0.0,
            absolute_error=max_error,
            relative_error=max_error,
            tolerance=0.20,  # 20% tolerance for rheology
            passed=max_error < 0.20,
            units="Pa·s",
            description="Non-Newtonian viscosity vs Casson model"
        )
        self.results.append(result_val)
        
        print(f"  Max viscosity error: {max_error*100:.2f}%")
        print(f"  Status: {'PASS' if result_val.passed else 'FAIL'}")
        
        return result_val
    
    def run_all_validations(self) -> Dict[str, bool]:
        """Run all 3D validation tests."""
        print("=" * 70)
        print("3D Flow Validation")
        print("=" * 70)
        
        results = {}
        results["poiseuille"] = self.validate_poiseuille_3d().passed
        results["bifurcation"] = self.validate_bifurcation_murray().passed
        results["wss"] = self.validate_wall_shear_stress().passed
        results["trifurcation"] = self.validate_trifurcation_flow().passed
        results["blood"] = self.validate_blood_rheology().passed
        
        self.generate_report()
        
        return results
    
    def generate_report(self) -> None:
        """Generate validation report."""
        print("\n" + "=" * 70)
        print("3D Validation Summary")
        print("=" * 70)
        
        passed = sum(1 for r in self.results if r.passed)
        total = len(self.results)
        
        print(f"\nTotal tests: {total}")
        print(f"Passed: {passed}")
        print(f"Failed: {total - passed}")
        print(f"Success rate: {passed/total*100:.1f}%" if total > 0 else "N/A")
        
        # Print detailed results
        print("\nDetailed Results:")
        print("-" * 70)
        for r in self.results:
            status = "PASS" if r.passed else "FAIL"
            print(f"\n{r.test_name}")
            print(f"  Status: {status}")
            print(f"  CFD-rs:   {r.cfrs_value:.6e} {r.units}")
            print(f"  Ref:      {r.reference_value:.6e} {r.units}")
            print(f"  Rel.err:  {r.relative_error*100:.4f}%")
        
        # Save JSON report
        report_file = self.output_dir / "validation_3d.json"
        with open(report_file, 'w') as f:
            json.dump([r.to_dict() for r in self.results], f, indent=2)
        print(f"\nReport saved to: {report_file}")
        
        # Generate plots if available
        if HAS_MATPLOTLIB:
            self._generate_plots()
    
    def _generate_plots(self) -> None:
        """Generate validation plots."""
        if not self.results:
            return
        
        fig, axes = plt.subplots(1, 2, figsize=(12, 5))
        
        # Error bar chart
        test_names = [r.test_name[:15] for r in self.results]
        rel_errors = [r.relative_error * 100 for r in self.results]
        tolerances = [r.tolerance * 100 for r in self.results]
        colors = ['green' if r.passed else 'red' for r in self.results]
        
        x = np.arange(len(test_names))
        axes[0].bar(x, rel_errors, color=colors, alpha=0.7)
        axes[0].plot(x, tolerances, 'b--', label='Tolerance', linewidth=2)
        axes[0].set_xlabel('Test Case')
        axes[0].set_ylabel('Relative Error (%)')
        axes[0].set_title('3D Validation Errors')
        axes[0].set_xticks(x)
        axes[0].set_xticklabels(test_names, rotation=45, ha='right')
        axes[0].legend()
        axes[0].set_yscale('log')
        
        # Pass/Fail pie chart
        passed_count = sum(1 for r in self.results if r.passed)
        failed_count = len(self.results) - passed_count
        axes[1].pie([passed_count, failed_count],
                   labels=['Passed', 'Failed'],
                   colors=['green', 'red'],
                   autopct='%1.1f%%')
        axes[1].set_title('3D Test Results')
        
        plt.tight_layout()
        plot_file = self.output_dir / "validation_3d.png"
        plt.savefig(plot_file, dpi=150, bbox_inches='tight')
        print(f"Plots saved to: {plot_file}")
        plt.close()


def main():
    parser = argparse.ArgumentParser(description="3D Flow Validation")
    parser.add_argument("--output", default="external_validation/results")
    parser.add_argument("--test",
                       choices=["all", "poiseuille", "bifurcation", "wss", "trifurcation", "blood"],
                       default="all")
    
    args = parser.parse_args()
    
    validator = Validator3D(output_dir=args.output)
    
    if args.test == "all":
        results = validator.run_all_validations()
    else:
        test_map = {
            "poiseuille": validator.validate_poiseuille_3d,
            "bifurcation": validator.validate_bifurcation_murray,
            "wss": validator.validate_wall_shear_stress,
            "trifurcation": validator.validate_trifurcation_flow,
            "blood": validator.validate_blood_rheology,
        }
        if args.test in test_map:
            result = test_map[args.test]()
            validator.generate_report()
            results = {args.test: result.passed}
    
    failed = sum(1 for passed in results.values() if not passed)
    return 1 if failed > 0 else 0


if __name__ == "__main__":
    sys.exit(main())
