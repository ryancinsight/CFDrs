# xtask - CFD-rs Build Automation

Automated build and validation system for CFD-rs Python bindings (cfd-python).

## Overview

The xtask system automates the complete workflow for building, testing, and validating CFD-rs:

1. **Build** cfd-python wheels with maturin
2. **Setup** Python virtual environment
3. **Install** dependencies (numpy, scipy, matplotlib, FEniCS)
4. **Validate** CFD implementations against analytical solutions and Python packages
5. **Report** results with plots and metrics

## Quick Start

### Complete Workflow (Recommended)

```bash
# Run everything: build, setup, install, validate
cargo xtask all --plot
```

This will:
- ✅ Check for compilation errors
- ✅ Build cfd-python wheel (release mode)
- ✅ Create Python virtual environment (.venv)
- ✅ Install all dependencies
- ✅ Run analytical validations
- ✅ Generate validation report with plots

### Step-by-Step Workflow

```bash
# 1. Check for compilation errors
cargo xtask check

# 2. Build cfd-python wheel
cargo xtask build-wheel --release

# 3. Setup virtual environment
cargo xtask setup-venv

# 4. Install dependencies and cfd-python
cargo xtask install-deps

# 5. Run validations
cargo xtask validate --plot --category analytical
```

## Commands

### `cargo xtask check`

Check if Rust code compiles without errors.

```bash
cargo xtask check
```

**What it does:**
- Runs `cargo build --workspace --all-features`
- Reports compilation errors if any
- Identifies missing trait bounds (common issue)

**Output:**
```
🔍 Checking for compilation errors...

✅ No compilation errors found!
Ready to build cfd-python wheels.
```

### `cargo xtask build-wheel`

Build cfd-python Python wheel with maturin.

```bash
# Release build (optimized)
cargo xtask build-wheel --release

# Debug build (faster, for development)
cargo xtask build-wheel --no-release
```

**What it does:**
- Checks compilation first
- Installs maturin if needed
- Builds wheel in `crates/cfd-python/target/wheels/`
- Lists available wheels

**Output:**
```
🔨 Building cfd-python wheel...

Building release wheel...
✅ Wheel built successfully!
Location: crates/cfd-python/target/wheels

Available wheels:
  - cfd-python-0.1.0-cp311-cp311-win_amd64.whl
```

### `cargo xtask setup-venv`

Create Python virtual environment.

```bash
# Default: use 'python' and create '.venv'
cargo xtask setup-venv

# Custom Python
cargo xtask setup-venv --python python3.11

# Custom path
cargo xtask setup-venv --path my-venv
```

**What it does:**
- Checks Python version
- Creates virtual environment
- Shows activation command

**Output:**
```
🐍 Setting up Python virtual environment...

Using: Python 3.11.5
Creating virtual environment...
✅ Virtual environment created at: .venv

Activate with:
  .venv\Scripts\activate
```

### `cargo xtask install-deps`

Install Python dependencies and cfd-python.

```bash
# Install all dependencies
cargo xtask install-deps

# Use custom venv path
cargo xtask install-deps --venv my-venv

# Flag to install FEniCS (shows instructions)
cargo xtask install-deps --with-fenics
```

**What it does:**
- Upgrades pip
- Installs maturin
- Installs requirements from `validation/requirements.txt`:
  - numpy, scipy, matplotlib, pandas
  - seaborn, plotly (visualization)
  - pytest (testing)
  - tabulate (reports)
- Installs cfd-python wheel if available
- Shows FEniCS installation instructions

**Output:**
```
📦 Installing Python dependencies...

Upgrading pip...
Installing maturin...
Installing validation dependencies...
Installing cfd-python...
✅ cfd-python installed!

✅ Dependencies installed!
```

### `cargo xtask install-fenics`

Install FEniCS for cross-validation.

```bash
# Install via conda (recommended)
cargo xtask install-fenics --use-conda

# Show pip instructions
cargo xtask install-fenics --no-use-conda
```

**What it does:**
- Checks if conda is available
- Creates conda environment `cfdrs-validation`
- Installs FEniCS from conda-forge
- Shows activation instructions

**Requirements:**
- Conda or Miniconda installed
- Internet connection (downloads ~500MB)

**Output:**
```
🧮 Installing FEniCS for validation...

Creating conda environment 'cfdrs-validation'...
This may take several minutes...

✅ FEniCS installed in conda environment!

Activate with:
  conda activate cfdrs-validation

Then install cfd-python:
  cargo xtask install-deps --venv $(conda info --base)/envs/cfdrs-validation
```

**Alternative (Docker):**
```bash
# Use pre-built FEniCS Docker image
docker run -ti quay.io/fenicsproject/stable:latest
```

### `cargo xtask validate`

Run validation suite.

```bash
# Run all validations with plots
cargo xtask validate --plot

# Quick mode (coarse grids, faster)
cargo xtask validate --quick

# Specific category
cargo xtask validate --category analytical
cargo xtask validate --category fenics
cargo xtask validate --category literature

# Custom venv
cargo xtask validate --venv my-venv
```

**What it does:**
- Checks if cfd-python is installed
- Runs validation scripts:
  - `analytical/poiseuille_2d.py` - 2D channel flow
  - `analytical/bernoulli_venturi.py` - Venturi validation
  - `fenics_comparison/poiseuille_fenics.py` - FEniCS comparison
- Generates report at `validation/reports/validation_summary.md`
- Creates plots in `validation/reports/figures/`

**Output:**
```
🧪 Running validation suite...

Running: python run_all_validations.py --plot --category analytical

================================================================================
ANALYTICAL SOLUTION VALIDATIONS
================================================================================

Running 2D Poiseuille validation...
✓ 2D Poiseuille completed

Running Venturi validation...
✓ Venturi completed

✅ Validation complete!
Report: validation/reports/validation_summary.md
```

### `cargo xtask all`

Complete workflow: build, setup, install, validate.

```bash
# Full workflow with plots
cargo xtask all --plot

# Include FEniCS installation
cargo xtask all --with-fenics --plot
```

**What it does:**
1. Check compilation
2. Build wheel (release)
3. Setup venv (if needed)
4. Install dependencies
5. Install FEniCS (optional)
6. Run validation suite
7. Generate report

**Output:**
```
🚀 Running complete validation workflow...

Step 1/6: Checking compilation...
Step 2/6: Building cfd-python wheel...
Step 3/6: Setting up virtual environment...
Step 4/6: Installing dependencies...
Step 5/6: Skipping FEniCS installation (use --with-fenics to enable)
Step 6/6: Running validation suite...

🎉 Complete workflow finished successfully!

Next steps:
  - Review report: validation/reports/validation_summary.md
  - View plots: validation/reports/figures/
  - Run full validation: cargo xtask validate --plot
```

### `cargo xtask clean`

Clean build artifacts and virtual environment.

```bash
cargo xtask clean
```

**What it does:**
- Runs `cargo clean`
- Removes `crates/cfd-python/target/wheels/`
- Removes `.venv/`
- Removes `validation/reports/`

**Output:**
```
🧹 Cleaning build artifacts...

Cleaning Rust build artifacts...
Removing cfd-python wheels...
Removing virtual environment...
Removing validation reports...

✅ Clean complete!
```

## Cargo Aliases

Convenient shortcuts in `.cargo/config.toml`:

```bash
# Same as: cargo xtask build-wheel --release
cargo build-python

# Same as: cargo xtask all --plot
cargo setup-validation

# Development
cargo dev              # Build workspace
cargo test-all         # Test all crates
cargo check-all        # Check all crates
cargo fmt-all          # Format code
cargo clippy-all       # Run clippy
```

## Validation Workflow

### 1. Analytical Validations (No External Dependencies)

Compares CFD-rs against exact mathematical solutions:

```bash
cargo xtask validate --category analytical --plot
```

**Tests:**
- **2D Poiseuille Flow:** Validates against `u(y) = -(dp/dx)/(2μ) · y(H-y)`
  - Acceptance: L2 error < 1%
- **Venturi Flow:** Validates against Bernoulli equation
  - Acceptance: Cp error < 10%

**Requirements:**
- Python 3.11+
- numpy, scipy, matplotlib
- cfd-python wheel

**Time:** ~2-5 minutes

### 2. FEniCS Comparison (Requires FEniCS)

Cross-validates against established FEM package:

```bash
# Install FEniCS first
cargo xtask install-fenics --use-conda
conda activate cfdrs-validation

# Run comparison
cargo xtask validate --category fenics --plot
```

**Tests:**
- **2D Poiseuille:** Compares CFD-rs vs FEniCS P2-P1 solution
  - Acceptance: Error < 5%
- **3D Bifurcation:** Wall shear stress comparison
  - Acceptance: WSS error < 10%

**Requirements:**
- FEniCS installed via conda
- cfd-python wheel
- ~1GB disk space

**Time:** ~10-20 minutes (FEniCS install) + ~5-10 minutes (validation)

### 3. Literature Benchmarks (Future)

Validates against published data:

```bash
cargo xtask validate --category literature --plot
```

**Planned tests:**
- Huo & Kassab (2012) bifurcation WSS
- ISO 5167 Venturi standards
- Stroock et al. (2002) serpentine mixer

## Troubleshooting

### Issue: Compilation errors

```
❌ Compilation errors detected!

Common issues:
  - Missing ToPrimitive trait bounds in cfd-1d/cfd-3d
  - Type conversion issues in generic implementations
```

**Solution:**
1. Check error details: `cargo build --workspace`
2. Add missing trait bounds:
   ```rust
   impl<T: RealField + Copy + ToPrimitive> MyStruct<T>
   ```
3. Run: `cargo xtask check` to verify

### Issue: Maturin not found

```
error: failed to run custom build command for `pyo3-ffi`
```

**Solution:**
```bash
pip install maturin
cargo xtask build-wheel
```

### Issue: Python version mismatch

```
error: the configured Python interpreter version (3.13) is newer than PyO3's maximum supported version (3.12)
```

**Solution:**
- Use Python 3.11 or 3.12
- Or set: `PYO3_USE_ABI3_FORWARD_COMPATIBILITY=1`

### Issue: cfd-python import fails

```python
ImportError: DLL load failed while importing cfd-python
```

**Solution:**
1. Check Python version matches build
2. Rebuild wheel: `cargo xtask build-wheel --release`
3. Reinstall: `cargo xtask install-deps`

### Issue: FEniCS not available

```
ModuleNotFoundError: No module named 'dolfin'
```

**Solution:**
```bash
# Install via conda
cargo xtask install-fenics --use-conda
conda activate cfdrs-validation

# Or use Docker
docker run -ti quay.io/fenicsproject/stable:latest
```

### Issue: Virtual environment exists

```
⚠️  Virtual environment already exists at: .venv
```

**Solution:**
```bash
# Remove and recreate
cargo xtask clean
cargo xtask setup-venv
```

## Directory Structure

```
CFDrs/
├── .cargo/
│   └── config.toml          # Cargo aliases
├── xtask/
│   ├── Cargo.toml           # xtask dependencies
│   ├── src/
│   │   └── main.rs          # Build automation
│   └── README.md            # This file
├── crates/cfd-python/
│   ├── src/                 # Python bindings
│   ├── target/wheels/       # Built wheels (generated)
│   └── Cargo.toml
├── validation/
│   ├── analytical/          # Analytical validations
│   ├── fenics_comparison/   # FEniCS comparisons
│   ├── requirements.txt     # Python dependencies
│   ├── run_all_validations.py
│   └── reports/             # Generated reports
│       ├── validation_summary.md
│       └── figures/         # Plots
└── .venv/                   # Virtual environment (generated)
```

## Development Workflow

### Typical development cycle:

```bash
# 1. Make changes to CFD-rs Rust code
vim crates/cfd-2d/src/solvers/simple.rs

# 2. Check compilation
cargo xtask check

# 3. Rebuild Python bindings
cargo xtask build-wheel --release

# 4. Reinstall in venv
cargo xtask install-deps

# 5. Run validations to verify
cargo xtask validate --quick

# 6. Full validation before commit
cargo xtask validate --plot
```

### Adding new validation tests:

```bash
# 1. Create validation script
vim validation/analytical/my_new_test.py

# 2. Add to run_all_validations.py
vim validation/run_all_validations.py

# 3. Run to test
.venv/Scripts/python validation/analytical/my_new_test.py

# 4. Add to suite
cargo xtask validate --category analytical
```

## CI/CD Integration

### GitHub Actions example:

```yaml
name: Validation

on: [push, pull_request]

jobs:
  validate:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v3
      
      - name: Setup Rust
        uses: actions-rs/toolchain@v1
        with:
          toolchain: stable
      
      - name: Setup Python
        uses: actions/setup-python@v4
        with:
          python-version: '3.11'
      
      - name: Run validation
        run: cargo xtask all --plot
      
      - name: Upload results
        uses: actions/upload-artifact@v3
        with:
          name: validation-report
          path: validation/reports/
```

## Environment Variables

### Optional configuration:

```bash
# Use specific Python version
export PYO3_PYTHON=/usr/bin/python3.11
cargo xtask build-wheel

# Skip Python version check
export PYO3_USE_ABI3_FORWARD_COMPATIBILITY=1
cargo xtask build-wheel

# Custom virtual environment path
cargo xtask setup-venv --path ~/my-custom-venv
```

## Performance

### Build times (approximate):

| Task | Time (Debug) | Time (Release) |
|------|--------------|----------------|
| Check compilation | ~30s | — |
| Build wheel | ~60s | ~3-5 min |
| Setup venv | ~10s | — |
| Install deps | ~30s | — |
| Install FEniCS | — | ~10-15 min |
| Run validation (quick) | ~1 min | — |
| Run validation (full) | ~3-5 min | — |
| **Complete workflow** | **~8-10 min** | **~15-25 min** |

### Disk space:

- Rust build artifacts: ~2-3 GB
- Python wheels: ~50-100 MB
- Virtual environment: ~200-500 MB
- FEniCS (conda): ~1-2 GB
- **Total:** ~3-6 GB

## Next Steps

1. **First time setup:**
   ```bash
   cargo xtask all --plot
   ```

2. **Daily development:**
   ```bash
   cargo xtask check
   cargo xtask build-wheel --release
   cargo xtask validate --quick
   ```

3. **Before committing:**
   ```bash
   cargo xtask validate --plot
   git add validation/reports/
   git commit -m "Update validation results"
   ```

4. **For publications:**
   ```bash
   cargo xtask all --with-fenics --plot
   # Review validation/reports/validation_summary.md
   ```

## References

- [xtask pattern](https://github.com/matklad/cargo-xtask)
- [maturin documentation](https://www.maturin.rs/)
- [FEniCS installation](https://fenicsproject.org/download/)
- [PyO3 guide](https://pyo3.rs/)

---

*xtask v0.1.0 - CFD-rs Build Automation*
