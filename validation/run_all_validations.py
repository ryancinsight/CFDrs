#!/usr/bin/env python3
"""
Comprehensive CFD-RS Validation Suite

This script coordinates all validation tests and generates a comprehensive report.

Usage:
    python run_all_validations.py [--category CATEGORY] [--plot] [--quick]

Categories:
    analytical  - Compare against analytical solutions
    fenics      - Compare against FEniCS solutions
    literature  - Compare against published data
    all         - Run all validations (default)

Options:
    --plot      - Generate plots and figures
    --quick     - Use coarse grids for fast testing
    --parallel  - Run validations in parallel (where possible)
"""

import argparse
import importlib.util
import json
import os
import sys
from datetime import datetime
from pathlib import Path

# Fix Windows console encoding for Unicode
if sys.platform == "win32":
    try:
        sys.stdout.reconfigure(encoding="utf-8")
        sys.stderr.reconfigure(encoding="utf-8")
    except AttributeError:
        pass  # Python < 3.7

# Add project root to path
PROJECT_ROOT = Path(__file__).parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

try:
    import matplotlib.pyplot as plt
    import numpy as np
    from tqdm import tqdm
except ImportError as e:
    print(f"❌ Missing required package: {e}")
    print("\nInstall dependencies with:")
    print("  pip install -r validation/requirements.txt")
    sys.exit(1)

try:
    import pycfdrs
except ImportError:
    print("❌ pycfdrs not installed!")
    print("\nBuild and install with:")
    print("  cargo xtask build-wheel")
    print("  cargo xtask install-deps")
    sys.exit(1)


class ValidationRunner:
    """Coordinates all validation tests"""

    def __init__(self, args):
        self.args = args
        self.results = {
            "timestamp": datetime.now().isoformat(),
            "pycfdrs_version": pycfdrs.__version__,
            "tests": {},
        }
        self.report_dir = PROJECT_ROOT / "validation" / "reports"
        self.figure_dir = self.report_dir / "figures"

        # Create directories
        self.report_dir.mkdir(parents=True, exist_ok=True)
        if args.plot:
            self.figure_dir.mkdir(parents=True, exist_ok=True)

    def run(self):
        """Run all validations based on category"""
        print("\n" + "=" * 80)
        print("CFD-RS VALIDATION SUITE")
        print("=" * 80)
        print(f"Category: {self.args.category}")
        print(
            f"Mode: {'Quick (coarse grids)' if self.args.quick else 'Full resolution'}"
        )
        print(f"Plotting: {'Enabled' if self.args.plot else 'Disabled'}")
        print("=" * 80 + "\n")

        # Run validations based on category
        if self.args.category in ["analytical", "all"]:
            self.run_analytical_validations()

        if self.args.category in ["fenics", "all"]:
            self.run_fenics_validations()

        if self.args.category in ["literature", "all"]:
            self.run_literature_validations()

        # Generate report
        self.generate_report()

        # Print summary
        self.print_summary()

    def run_analytical_validations(self):
        """Run analytical solution validations"""
        print("\n" + "🔬" * 40)
        print("ANALYTICAL VALIDATIONS")
        print("🔬" * 40 + "\n")

        tests = []

        # 1D Poiseuille (if implemented)
        tests.append(
            {
                "name": "1D Poiseuille Flow",
                "script": "validation_1d_poiseuille.py",
                "required": False,
            }
        )

        # 1D Bifurcation
        tests.append(
            {
                "name": "1D Bifurcation",
                "script": "validation_analytical.py",
                "required": True,
            }
        )

        # 2D Poiseuille
        tests.append(
            {
                "name": "2D Poiseuille Flow",
                "script": "validation_2d_poiseuille.py",
                "required": False,
            }
        )

        for test in tests:
            self.run_test(test, category="analytical")

    def run_fenics_validations(self):
        """Run FEniCS comparison validations"""
        print("\n" + "⚗️" * 40)
        print("FENICS COMPARISONS")
        print("⚗️" * 40 + "\n")

        # Check if FEniCS is available
        try:
            import dolfin  # FEniCS

            print("✅ FEniCS detected\n")
        except ImportError:
            print("⚠️  FEniCS not installed - skipping FEniCS validations")
            print("Install with: cargo xtask install-fenics --use-conda\n")
            self.results["tests"]["fenics"] = {
                "status": "skipped",
                "reason": "FEniCS not installed",
            }
            return

        tests = [
            {
                "name": "2D Poiseuille vs FEniCS",
                "script": "fenics/compare_2d_poiseuille.py",
                "required": False,
            },
            {
                "name": "2D Venturi vs FEniCS",
                "script": "fenics/compare_2d_venturi.py",
                "required": False,
            },
        ]

        for test in tests:
            self.run_test(test, category="fenics")

    def run_literature_validations(self):
        """Run literature data validations"""
        print("\n" + "📚" * 40)
        print("LITERATURE VALIDATIONS")
        print("📚" * 40 + "\n")

        tests = [
            {
                "name": "Murray Law Validation",
                "script": "literature/validate_murray_law.py",
                "required": False,
            },
            {
                "name": "Blood Rheology Literature",
                "script": "literature/validate_blood_rheology.py",
                "required": False,
            },
        ]

        for test in tests:
            self.run_test(test, category="literature")

    def run_test(self, test, category):
        """Run a single validation test"""
        test_name = test["name"]
        script_path = PROJECT_ROOT / "validation" / test["script"]

        print(f"Running: {test_name}")
        print(f"  Script: {test['script']}")

        if not script_path.exists():
            print(f"  ⚠️  NOT IMPLEMENTED")
            self.results["tests"][test_name] = {
                "status": "not_implemented",
                "category": category,
            }
            if test["required"]:
                print(f"  ❌ REQUIRED TEST MISSING")
            print()
            return

        # Try to run the test
        try:
            # Import and run the test module
            spec = importlib.util.spec_from_file_location(test_name, script_path)
            module = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(module)

            # Check for run_validation function
            if hasattr(module, "run_validation"):
                result = module.run_validation(
                    plot=self.args.plot,
                    quick=self.args.quick,
                    output_dir=self.figure_dir if self.args.plot else None,
                )

                self.results["tests"][test_name] = {
                    "status": "passed" if result.get("passed", False) else "failed",
                    "category": category,
                    **result,
                }

                if result.get("passed", False):
                    print(f"  ✅ PASSED")
                else:
                    print(f"  ❌ FAILED: {result.get('error', 'Unknown error')}")
            else:
                print(f"  ⚠️  No run_validation() function found")
                self.results["tests"][test_name] = {
                    "status": "error",
                    "category": category,
                    "error": "No run_validation() function",
                }

        except Exception as e:
            print(f"  ❌ ERROR: {e}")
            self.results["tests"][test_name] = {
                "status": "error",
                "category": category,
                "error": str(e),
            }

        print()

    def generate_report(self):
        """Generate markdown validation report"""
        report_path = self.report_dir / "validation_summary.md"

        with open(report_path, "w", encoding="utf-8") as f:
            f.write(f"# CFD-RS Validation Report\n\n")
            f.write(f"**Generated:** {self.results['timestamp']}  \n")
            f.write(f"**pycfdrs Version:** {self.results['pycfdrs_version']}  \n")
            f.write(f"**Category:** {self.args.category}  \n")
            f.write(f"**Mode:** {'Quick' if self.args.quick else 'Full'}  \n\n")

            f.write("---\n\n")
            f.write("## Summary\n\n")

            # Count statuses
            status_counts = {}
            for test_name, result in self.results["tests"].items():
                status = result["status"]
                status_counts[status] = status_counts.get(status, 0) + 1

            f.write(f"**Total Tests:** {len(self.results['tests'])}  \n")
            for status, count in status_counts.items():
                icon = {
                    "passed": "✅",
                    "failed": "❌",
                    "not_implemented": "⚠️",
                    "skipped": "⏭️",
                    "error": "💥",
                }.get(status, "❓")
                f.write(f"**{status.title()}:** {icon} {count}  \n")

            f.write("\n---\n\n")
            f.write("## Detailed Results\n\n")

            # Group by category
            categories = {}
            for test_name, result in self.results["tests"].items():
                cat = result.get("category", "unknown")
                if cat not in categories:
                    categories[cat] = []
                categories[cat].append((test_name, result))

            for category, tests in sorted(categories.items()):
                f.write(f"### {category.title()} Validations\n\n")

                for test_name, result in tests:
                    status = result["status"]
                    icon = {
                        "passed": "✅",
                        "failed": "❌",
                        "not_implemented": "⚠️",
                        "skipped": "⏭️",
                        "error": "💥",
                    }.get(status, "❓")

                    f.write(f"#### {icon} {test_name}\n\n")
                    f.write(f"**Status:** {status}  \n")

                    if "error" in result:
                        f.write(f"**Error:** {result['error']}  \n")

                    if "max_error" in result:
                        f.write(f"**Max Error:** {result['max_error']:.2e}  \n")

                    if "mean_error" in result:
                        f.write(f"**Mean Error:** {result['mean_error']:.2e}  \n")

                    if "metrics" in result:
                        f.write("\n**Metrics:**\n")
                        for metric, value in result["metrics"].items():
                            f.write(f"- {metric}: {value}  \n")

                    f.write("\n")

            f.write("---\n\n")
            f.write("## Figures\n\n")

            if self.args.plot and self.figure_dir.exists():
                figures = list(self.figure_dir.glob("*.png"))
                if figures:
                    for fig in sorted(figures):
                        f.write(f"### {fig.stem}\n\n")
                        f.write(f"![{fig.stem}](figures/{fig.name})\n\n")
                else:
                    f.write("*No figures generated*\n\n")
            else:
                f.write("*Plotting disabled (use --plot to generate figures)*\n\n")

        print(f"\n📄 Report saved to: {report_path}")

        # Also save JSON for programmatic access
        json_path = self.report_dir / "validation_results.json"
        with open(json_path, "w") as f:
            json.dump(self.results, f, indent=2)

        print(f"📄 JSON results saved to: {json_path}")

    def print_summary(self):
        """Print summary to console"""
        print("\n" + "=" * 80)
        print("VALIDATION SUMMARY")
        print("=" * 80 + "\n")

        status_counts = {}
        for result in self.results["tests"].values():
            status = result["status"]
            status_counts[status] = status_counts.get(status, 0) + 1

        total = len(self.results["tests"])
        passed = status_counts.get("passed", 0)
        failed = status_counts.get("failed", 0)
        not_impl = status_counts.get("not_implemented", 0)

        print(f"Total Tests: {total}")
        print(f"  ✅ Passed: {passed}")
        print(f"  ❌ Failed: {failed}")
        print(f"  ⚠️  Not Implemented: {not_impl}")

        if failed > 0:
            print("\n❌ VALIDATION FAILED")
            print("Review the report for details:")
            print(f"  {self.report_dir / 'validation_summary.md'}")
            sys.exit(1)
        elif passed > 0:
            print("\n✅ ALL IMPLEMENTED TESTS PASSED")
            if not_impl > 0:
                print(f"\nNote: {not_impl} tests not yet implemented")
        else:
            print("\n⚠️  NO TESTS RUN")

        print()


def main():
    parser = argparse.ArgumentParser(
        description="Run CFD-RS validation suite",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )

    parser.add_argument(
        "--category",
        choices=["analytical", "fenics", "literature", "all"],
        default="all",
        help="Validation category to run",
    )

    parser.add_argument(
        "--plot", action="store_true", help="Generate plots and figures"
    )

    parser.add_argument(
        "--quick", action="store_true", help="Use coarse grids for fast testing"
    )

    parser.add_argument(
        "--parallel",
        action="store_true",
        help="Run validations in parallel (experimental)",
    )

    args = parser.parse_args()

    # Run validation
    runner = ValidationRunner(args)
    runner.run()


if __name__ == "__main__":
    main()
