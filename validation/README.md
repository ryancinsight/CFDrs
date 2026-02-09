# Validation Suite

This directory contains validation scripts for the 3D CFD solvers.

## Prerequisites

- `pycfdrs` installed (`maturin develop` in project root)
- `numpy`
- `scipy` (optional, for analytical comparisons)

## Running validation

Run individual scripts:

```bash
python validation/validate_bifurcation.py
python validation/validate_trifurcation.py
python validation/validate_venturi.py
python validation/validate_serpentine.py
```

## Metrics

Each script checks:
- Geometric fidelity
- Mass conservation
- Physics consistency (pressure drops, velocity scaling)
- Solver convergence (via bindings)
