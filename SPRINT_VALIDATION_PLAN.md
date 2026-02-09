# Sprint 1.90.0: Cross-Package CFD Validation

**Objective**: Complete rigorous validation by comparing cfd-rs against established Python CFD packages

## Target External Packages
1. **DrZGan/Python_CFD** - Classic CFD teaching codes
2. **pmocz/cfd-comparison-python** - CFD comparison benchmarks  
3. **fluidsim** - Spectral fluid simulation framework

## Test Case Priority Matrix

| Test Case | Status | Target Package | Literature Reference |
|-----------|--------|----------------|---------------------|
| Poiseuille 2D | ✅ Working | Analytical | Hagen-Poiseuille |
| Lid-Driven Cavity | ❌ Skipped | Ghia et al. 1982 | Ghia et al. (1982) |
| Bifurcation Blood Flow | 📋 Ready | Murray's Law | Murray (1926) |
| Trifurcation Flow | 📋 Ready | Murray's Law | Murray (1926) |
| Serpentine Mixing | 📋 Ready | Dean Flow Literature | Dean (1927) |
| Venturi Throat | 📋 Ready | Analytical | Bernoulli |

## Implementation Tasks

### Phase 1: External Package Setup
- [ ] Clone/install DrZGan/Python_CFD
- [ ] Clone/install pmocz/cfd-comparison-python
- [ ] Verify fluidsim installation
- [ ] Create validation wrapper scripts

### Phase 2: Lid-Driven Cavity Validation (CRITICAL)
- [ ] Diagnose SIMPLEC convergence issue
- [ ] Fix solver convergence
- [ ] Validate against Ghia et al. (1982)
- [ ] Compare with Python_CFD implementations

### Phase 3: Blood-Specific Validation
- [ ] Implement Carreau-Yasuda rheology validation
- [ ] Validate bifurcation flow vs Murray's Law
- [ ] Validate trifurcation resistance ratios
- [ ] Compare with published blood flow data

### Phase 4: Geometric Validation
- [ ] Serpentine Dean number correlation
- [ ] Venturi throat pressure recovery
- [ ] Bifurcation flow division ratios

## Success Criteria
- All test cases pass with <5% error vs literature/analytical
- Cross-package agreement within 2%
- Zero placeholders, stubs, or simplified paths
