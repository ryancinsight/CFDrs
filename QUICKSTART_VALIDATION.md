# CFD-rs Validation Quick Start Guide

Get your CFD-rs implementation validated in 5 minutes or less.

## Prerequisites

- ✅ Rust toolchain (1.70+)
- ✅ Python 3.11 or 3.12
- ✅ Git
- ⚠️ ~5 GB disk space
- ⚠️ Internet connection

## Option 1: Complete Automated Workflow (Recommended)

Run everything with a single command:

```bash
cargo xtask all --plot
```

This will:
1. ✅ Check for compilation errors
2. ✅ Build pycfdrs Python package
3. ✅ Setup virtual environment
4. ✅ Install all dependencies
5. ✅ Run analytical validations
6. ✅ Generate validation report with plots

**Time:** ~10-15 minutes (first run)

**Output:**
- Report: `validation/reports/validation_summary.md`
- Plots: `validation/reports/figures/*.png`
- Virtual environment: `.venv/`

### View Results

```bash
# On Windows
notepad validation\reports\validation_summary.md

# On Linux/Mac
cat validation/reports/validation_summary.md

# View plots
explorer validation\reports\figures       # Windows
open validation/reports/figures           # Mac
xdg-open validation/reports/figures       # Linux
```

## Option 2: Step-by-Step (For Learning)

### Step 1: Check Compilation

```bash
cargo xtask check
```

**Expected output:**
```
🔍 Checking for compilation errors...
✅ No compilation errors found!
```

**If errors occur:**
- These are pre-existing issues in cfd-1d/cfd-3d
- Fix by adding `ToPrimitive` trait bounds
- See `VALIDATION_IMPLEMENTATION_GUIDE.md` for details

### Step 2: Build Python Package

```bash
cargo xtask build-wheel --release
```

**Expected output:**
```
🔨 Building pycfdrs wheel...
✅ Wheel built successfully!
Location: crates/pycfdrs/target/wheels

Available wheels:
  - pycfdrs-0.1.0-cp311-cp311-win_amd64.whl
```

**Time:** ~3-5 minutes

### Step 3: Setup Environment

```bash
cargo xtask setup-venv
```

**Expected output:**
```
🐍 Setting up Python virtual environment...
✅ Virtual environment created at: .venv

Activate with:
  .venv\Scripts\activate    (Windows)
  source .venv/bin/activate (Linux/Mac)
```

**Time:** ~10 seconds

### Step 4: Install Dependencies

```bash
cargo xtask install-deps
```

**Expected output:**
```
📦 Installing Python dependencies...
Installing pycfdrs...
✅ pycfdrs installed!
✅ Dependencies installed!
```

**Time:** ~30-60 seconds

### Step 5: Run Validation

```bash
cargo xtask validate --plot --category analytical
```

**Expected output:**
```
🧪 Running validation suite...

================================================================================
2D POISEUILLE FLOW VALIDATION
================================================================================
✓ L2 error < 1%
✓ L∞ error < 5%
✓ ALL VALIDATION CHECKS PASSED

================================================================================
VENTURI FLOW VALIDATION
================================================================================
✓ Cp error < 10%
✓ Pressure recovery > 80%
✓ ALL VALIDATION CHECKS PASSED

✅ Validation complete!
```

**Time:** ~2-5 minutes

## Option 3: Quick Test (Fastest)

Minimal validation for quick feedback:

```bash
# Quick validation with coarse grids
cargo xtask validate --quick --category analytical
```

**Time:** ~1-2 minutes

## Understanding the Results

### Validation Report

Location: `validation/reports/validation_summary.md`

**Key sections:**
```markdown
## Executive Summary

**Total Tests:** 2
**Passed:** 2 ✓
**Failed:** 0 ✗
**Success Rate:** 100%

## Test Results

| Category   | Test          | Status  | Key Metrics |
|------------|---------------|---------|-------------|
| Analytical | Poiseuille 2D | ✓ PASS  | L2: 0.45%   |
| Analytical | Venturi       | ✓ PASS  | Cp: 8.2%    |
```

### Validation Plots

Location: `validation/reports/figures/`

**Generated plots:**
- `poiseuille_2d_validation.png` - Velocity profile comparison
- `venturi_validation.png` - Pressure distribution and geometry

**What to look for:**
- ✅ CFD-rs line matches analytical line closely
- ✅ Error < 5% across the domain
- ✅ No unphysical oscillations

### Acceptance Criteria

| Metric | Threshold | Your Result |
|--------|-----------|-------------|
| L2 relative error | < 1% | Check report |
| L∞ relative error | < 5% | Check report |
| Mass conservation | < 10⁻⁶ | Check report |
| Cp error (Venturi) | < 10% | Check report |

## Advanced: FEniCS Cross-Validation

Compare against established FEM package for extra confidence.

### Install FEniCS

```bash
cargo xtask install-fenics --use-conda
```

**Time:** ~10-15 minutes (one-time)

### Run FEniCS Comparison

```bash
# Activate FEniCS environment
conda activate cfdrs-validation

# Install pycfdrs in conda env
pip install crates/pycfdrs/target/wheels/*.whl

# Run comparison
python validation/fenics_comparison/poiseuille_fenics.py
```

**Expected output:**
```
FEniCS vs CFD-rs: 2D POISEUILLE FLOW COMPARISON
FEniCS max velocity: 1.234e-03 m/s
CFD-rs max velocity: 1.230e-03 m/s
Relative error: 0.32%
✓ FEniCS accuracy < 1%
```

## Troubleshooting

### Problem: Compilation errors

**Error message:**
```
❌ Compilation errors detected!
  - Missing ToPrimitive trait bounds in cfd-1d/cfd-3d
```

**Solution:**
1. Check specific errors: `cargo build --workspace`
2. Fix trait bounds in affected files
3. Verify: `cargo xtask check`

**Workaround:** Skip to validation framework testing with mock data (validation scripts will run with analytical solutions only)

### Problem: maturin not found

**Error message:**
```
error: maturin not found
```

**Solution:**
```bash
pip install maturin
cargo xtask build-wheel
```

### Problem: Python version too new

**Error message:**
```
error: Python 3.13 is newer than PyO3's maximum supported version (3.12)
```

**Solution:**
- Use Python 3.11 or 3.12
- Or run: `setup-venv --python python3.11`

### Problem: Virtual environment already exists

**Error message:**
```
⚠️ Virtual environment already exists at: .venv
```

**Solution:**
```bash
# Clean and recreate
cargo xtask clean
cargo xtask all --plot
```

### Problem: pycfdrs import fails in Python

**Error in Python:**
```python
>>> import pycfdrs
ImportError: DLL load failed
```

**Solution:**
1. Verify Python version: `python --version`
2. Rebuild wheel: `cargo xtask build-wheel --release`
3. Reinstall: `cargo xtask install-deps`

## Common Use Cases

### Use Case 1: Daily Development

```bash
# Make changes to Rust code
vim crates/cfd-2d/src/solvers/simple.rs

# Quick validation
cargo xtask check
cargo xtask build-wheel --release
cargo xtask validate --quick
```

### Use Case 2: Before Committing

```bash
# Full validation with plots
cargo xtask validate --plot

# Review results
cat validation/reports/validation_summary.md

# Commit with results
git add validation/reports/
git commit -m "Validated CFD implementation"
```

### Use Case 3: Publication/Presentation

```bash
# Complete validation including FEniCS
cargo xtask all --with-fenics --plot

# Generate high-quality plots
python validation/analytical/poiseuille_2d.py

# Include in paper
# - validation/reports/validation_summary.md (results table)
# - validation/reports/figures/*.png (plots)
```

### Use Case 4: Continuous Integration

```bash
# In CI/CD pipeline
cargo xtask all --plot --no-with-fenics

# Check exit code
if [ $? -eq 0 ]; then
  echo "Validation passed!"
else
  echo "Validation failed!"
  exit 1
fi
```

## What Gets Validated

### 1. Analytical Solutions (Always)

**2D Poiseuille Flow:**
- Exact solution: `u(y) = -(dp/dx)/(2μ) · y(H-y)`
- Tests: Velocity profile, flow rate, wall shear stress
- Physics: Laminar flow between parallel plates

**Venturi Flow:**
- Exact solution: Bernoulli equation
- Tests: Pressure coefficient, pressure recovery
- Physics: Flow through converging-diverging channel

### 2. FEniCS Comparison (Optional)

**2D Poiseuille with FEM:**
- FEniCS P2-P1 Taylor-Hood elements
- Compares velocity fields point-by-point
- Tests consistency between FVM (CFD-rs) and FEM (FEniCS)

### 3. Blood Flow Models (Built-in)

**Rheology validation:**
- Casson model (yield stress)
- Carreau-Yasuda model (shear-thinning)
- Literature-validated parameters

## Performance Expectations

### Expected Results

Based on existing CFD-rs quality:

| Component | Expected Pass Rate | Confidence |
|-----------|-------------------|------------|
| 2D Poiseuille | 95-100% | HIGH |
| Venturi | 85-95% | MEDIUM-HIGH |
| FEniCS comparison | 90-100% | HIGH |
| **Overall** | **85-95%** | **HIGH** |

### Timing

| Task | First Run | Subsequent Runs |
|------|-----------|-----------------|
| Build wheel | 3-5 min | 2-3 min |
| Setup venv | 10 sec | — |
| Install deps | 30-60 sec | 10 sec |
| Run validation | 2-5 min | 2-5 min |
| **Total** | **~10-15 min** | **~5-10 min** |

## Next Steps

### Immediate

1. Run: `cargo xtask all --plot`
2. Review: `validation/reports/validation_summary.md`
3. Check plots in: `validation/reports/figures/`

### Short Term

1. Install FEniCS for cross-validation
2. Run all validation categories
3. Add custom validation tests

### Long Term

1. Integrate into CI/CD
2. Automate validation on every commit
3. Track validation metrics over time
4. Add literature benchmarks

## Getting Help

### Documentation

- **This file:** Quick start guide
- **`xtask/README.md`:** Detailed xtask commands
- **`VALIDATION_IMPLEMENTATION_GUIDE.md`:** Complete methodology
- **`validation/README.md`:** Validation framework overview

### Common Commands

```bash
# List all xtask commands
cargo xtask --help

# Get help for specific command
cargo xtask validate --help

# Check what's installed
.venv/Scripts/python -c "import pycfdrs; print(pycfdrs.__version__)"

# Run specific validation
.venv/Scripts/python validation/analytical/poiseuille_2d.py
```

### Debugging

```bash
# Verbose output
RUST_LOG=debug cargo xtask build-wheel

# Check virtual environment
.venv/Scripts/python -c "import sys; print(sys.executable)"

# Check installed packages
.venv/Scripts/pip list

# Manual validation run
.venv/Scripts/activate
cd validation
python run_all_validations.py --plot --category analytical
```

## Success Criteria

You'll know validation is working when:

✅ `cargo xtask all --plot` completes without errors  
✅ Report shows "All tests passed"  
✅ Plots show CFD-rs matching analytical solutions  
✅ L2 error < 1% for Poiseuille flow  
✅ Cp error < 10% for Venturi flow  

## Example Session

Complete workflow from start to finish:

```bash
# Clone repository (if needed)
git clone <repo-url>
cd CFDrs

# Run complete validation
cargo xtask all --plot

# Check results
cat validation/reports/validation_summary.md

# Expected output:
# ✓ ALL VALIDATION CHECKS PASSED
# Success Rate: 100%

# View plots
explorer validation\reports\figures

# Optional: Run with FEniCS
cargo xtask install-fenics --use-conda
conda activate cfdrs-validation
python validation/fenics_comparison/poiseuille_fenics.py

# Clean up
cargo xtask clean
```

---

**🎉 Ready to validate CFD-rs!**

Run: `cargo xtask all --plot`

---

*CFD-rs Validation Quick Start v1.0*  
*For detailed documentation, see `xtask/README.md`*
