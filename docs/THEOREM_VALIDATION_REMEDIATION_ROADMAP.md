# CFDrs Mathematical Rigor Remediation Roadmap

**Version**: 1.0.0  
**Date**: 2025-11-03  
**Based On**: Theorem Validation & Test Rigor Audit v1.0.0  
**Scope**: Sprint planning for mathematical completeness and test validation enhancement

---

## Executive Summary

This roadmap provides a systematic plan to elevate CFDrs from its current **mixed mathematical rigor (50-100%)** to **>85% across all components** within 12 sprints (~36 hours). The plan prioritizes **critical production blockers** (SIMPLE/SIMPLEC accuracy, hardcoded bugs) before addressing **theoretical documentation gaps** (linear solver convergence proofs, AMG theory).

### Remediation Metrics

**Current State** (from audit):
- Components with >90% rigor: 7/17 (41%)
- Critical issues (P0): 2
- High-priority gaps (P1): 3
- Production-ready components: 7/17 (41%)

**Target State** (Sprint 1.106.0):
- Components with >85% rigor: 17/17 (100%)
- Critical issues: 0
- High-priority gaps: 0
- Production-ready components: 17/17 (100%)

---

## Phase 1: Critical Production Blockers (Sprints 1.95.0-1.96.0, 8h)

### Sprint 1.95.0 - Critical Bug Fixes (4h)

**Objective**: Resolve showstopper bugs preventing model functionality

#### Task 1.1: Fix Rectangular Channel Hardcoded Reynolds (2h)

**Current Bug** (`crates/cfd-1d/src/resistance/models/rectangular.rs:108-109`):
```rust
let reynolds = T::from_f64(100.0).unwrap_or_else(|| T::one()); // HARDCODED!
```

**Fix**:
```rust
// Accept Reynolds from FlowConditions
let reynolds = match &self.conditions {
    Some(cond) => cond.reynolds_number(),
    None => return Err(Error::MissingFlowConditions),
};
```

**Validation**:
- Add 5 tests: Re ∈ {50, 100, 500, 1000, 2000}
- Verify Po scaling: square channel Po=56.91 at all Re
- Test aspect ratio dependency: α ∈ {0.2, 0.5, 1.0, 2.0, 5.0}

**Acceptance Criteria**:
- ✅ Model functional for arbitrary Reynolds numbers
- ✅ Tests pass for Re range [10, 2300]
- ✅ Aspect ratio tests show correct Po variation

---

#### Task 1.2: Document Rectangular Channel Formula Correction (2h)

**Current Documentation Issue**:
- Documented formula: R = (Po μ L) / (Re ρ A Dh)
- Actual implementation: R = f_fanning · L · ρ / (2 · A · Dh²)
- **These are not equivalent**

**Action**:
1. Derive implementation formula from Darcy-Weisbach:
   ```
   ΔP = f · (L/Dh) · (ρV²/2)
   R = ΔP / Q̇ = [f · (L/Dh) · (ρV²/2)] / (V·A) = f · ρ · L / (2·A·Dh)
   ```
   But code uses Dh² (line 111) - investigate if bug or different definition

2. Update documentation to match implementation OR fix implementation
3. Add derivation in comments showing Poiseuille number connection

**Acceptance Criteria**:
- ✅ Documentation matches implementation
- ✅ Derivation from first principles included
- ✅ Shah-London polynomial connection explained

---

### Sprint 1.96.0 - SIMPLE/SIMPLEC Accuracy Investigation (4h)

**Objective**: Root cause analysis for 30-35% Ghia cavity errors vs <8% literature standard

#### Task 1.3: Systematic SIMPLE Algorithm Debugging (4h)

**Current Failure Modes**:
- SIMPLEC Re=100: 30% L2 error (target <5-8%)
- SIMPLEC Re=400: 35% L2 error (target <8%)
- PIMPLE Re=100: 35% L2 error (target <5%)

**Debugging Protocol**:

1. **Pressure Correction Validation** (1h)
   - Verify Poisson equation discretization: ∇²p' = (ρ/Δt)∇·u*
   - Check boundary conditions: proper Neumann BC implementation?
   - Test reference pressure point: Does (1,1) hardcoding affect solution?

2. **Momentum Discretization Audit** (1h)
   - Verify convection scheme: upwind vs central differencing
   - Check diffusion discretization: proper second derivatives?
   - Validate source term assembly: pressure gradient included?

3. **Rhie-Chow Effectiveness** (1h)
   - Compare solutions with/without Rhie-Chow interpolation
   - Check face velocity formula: uf = ūf + df·[(∇p)P - (∇p)f]
   - Validate momentum coefficient interpolation: arithmetic vs harmonic mean?

4. **Under-Relaxation Tuning** (1h)
   - Test αu ∈ {0.5, 0.7, 0.9}, αp ∈ {0.1, 0.3, 0.5}
   - Monitor convergence behavior: monotonic residual decrease?
   - Compare with Patankar (1980) recommended values for lid-driven cavity

**Deliverable**:
- Technical report identifying root cause(s)
- Proposed fixes with theoretical justification
- If architectural issue: alternative algorithm recommendation

**Acceptance Criteria**:
- ✅ Root cause(s) documented with evidence
- ✅ Fix path identified OR decision to deprecate algorithm
- ✅ Stakeholder approval for next sprint implementation

---

## Phase 2: High-Priority Theoretical Gaps (Sprints 1.97.0-1.99.0, 12h)

### Sprint 1.97.0 - Entrance Effects Model Remediation (6h)

**Objective**: Provide theoretical foundation and validation for entrance loss model

#### Task 2.1: Derive Entrance Loss Formulas (2h)

**Current Gap**: Empirical correlations stated without derivation

**Action**:

1. **Sudden Contraction Derivation** (1h)
   - Start from Bernoulli + momentum balance
   - Show how K = (1 - A₁/A₂)² emerges from streamline contraction
   - Explain Reynolds correction: (1 + C/Re) from viscous losses
   - Justify C ∈ [0.1, 0.2] range with experimental data

2. **Developing Flow Length Theory** (1h)
   - Derive laminar: Lh/Dh ≈ 0.05Re from boundary layer growth
   - Cite Shah & London (1978) exact solution
   - Explain turbulent: Lh/Dh ≈ 4.4 Re^(1/6) from empirical correlation
   - Reference White (2006) §6.11 experimental data

**Deliverable**:
- Complete mathematical derivations in comments (lines 11-38)
- Literature citations with specific equation numbers
- Validity range documentation

---

#### Task 2.2: Justify or Replace Hardcoded Coefficients (2h)

**Current Issues**:
- Laminar entrance K=0.01 (line 172) - no source cited
- Smooth contraction regime split at Re=10⁴ - arbitrary?

**Action**:

1. **Literature Search** (1h)
   - Idelchik (1994) Handbook §§ 4.1-4.3 for entrance coefficients
   - White (2006) Table 6.4 for developing flow data
   - Find authoritative values for:
     - Laminar smooth contraction K(Re, A₁/A₂)
     - Turbulent threshold Re_transition

2. **Implementation Update** (1h)
   - Replace K=0.01 with literature value OR derive from theory
   - Document Re=10⁴ threshold justification OR implement smooth transition
   - Add validity checks: warn if conditions outside validated range

**Deliverable**:
- All magic numbers removed or justified
- Literature references for all coefficients
- Validity range enforcement

---

#### Task 2.3: Implement Entrance Effects Test Suite (2h)

**Current Gap**: Zero test coverage

**Test Cases**:

1. **Sudden Contraction Validation** (30min)
   - Test matrix: A₁/A₂ ∈ {0.2, 0.5, 0.8}, Re ∈ {100, 1000, 10000}
   - Compare K against Idelchik (1994) Table 4-1
   - Tolerance: ±10% (empirical correlation uncertainty)

2. **Smooth Contraction Validation** (30min)
   - Test: A₁/A₂ ∈ {0.3, 0.6, 0.9}
   - Verify K = 0.05 + 0.19·(A₁/A₂) formula
   - Compare with experimental data if available

3. **Developing Flow Length** (30min)
   - Laminar: Lh/Dh vs 0.05·Re for Re ∈ {100, 500, 1000}
   - Turbulent: Lh/Dh vs 4.4·Re^(1/6) for Re ∈ {5000, 10000, 50000}
   - Tolerance: ±20% (correlation scatter)

4. **Dimensional Analysis** (30min)
   - Verify R_entrance has units [Pa·s·m³]
   - Check conversion from K to R is dimensionally consistent
   - Test: changing units doesn't affect dimensionless results

**Acceptance Criteria**:
- ✅ 4 test categories implemented
- ✅ All tests passing with documented tolerances
- ✅ Experimental data cited for validation

---

### Sprint 1.98.0 - DES Shielding Function Enhancement (4h)

**Objective**: Implement proper Spalart et al. (2006) DDES formulation

#### Task 2.4: Implement Full DDES Shielding Function (2h)

**Current Simplified Form** (`crates/cfd-2d/src/physics/turbulence/des.rs:182-206`):
```rust
fd = 1 - tanh[(-dw/(0.5·Δ))]  // Incorrect
```

**Spalart et al. (2006) Correct Form**:
```rust
// Equation 18 from Spalart et al. (2006) AIAA Journal 44(11):2431-2441
let r_d = (nu_t + nu) / (kappa.powi(2) * dw.powi(2) * sqrt(0.5 * (S.powi(2) + Omega.powi(2))));
let f_d = 1.0 - tanh((8.0 * r_d).powi(3));

// CDES typically 0.65 for k-ω, 0.61 for Spalart-Allmaras
let l_des = if f_d < 1.0 {
    min(l_rans, (1.0 - f_d) * l_rans + f_d * cdes * delta)
} else {
    l_rans
};
```

**Implementation Steps**:
1. Compute strain rate magnitude S and vorticity Ω (reuse from S-A model)
2. Calculate wall distance dw efficiently (consider Eikonal equation solver)
3. Implement rd = (νt+ν)/(κ²dw²·√(S²+Ω²))
4. Apply shielding: fd = 1 - tanh[(8·rd)³]
5. Blend length scales: LDES = (1-fd)·LRANS + fd·CDES·Δ

**Deliverable**:
- Complete DDES shielding function implementation
- Comments citing Spalart et al. (2006) Eq. 18
- Unit tests for fd ∈ [0,1] in boundary layer and freestream

---

#### Task 2.5: Validate DES Against Literature Benchmark (2h)

**Benchmark**: Flat plate boundary layer with forced separation

**Setup**:
- Domain: Flat plate with adverse pressure gradient
- Reynolds: Re_θ ≈ 10,000 (transition to separation)
- Mesh: y+<1 in viscous sublayer, Δx+≈50-100 in LES region

**Validation Metrics**:
1. **Shielding function profile**:
   - fd ≈ 0 in boundary layer (RANS mode)
   - fd ≈ 1 in separated region (LES mode)
   - Smooth transition zone

2. **Separation point**:
   - Compare with DNS/experimental data
   - Tolerance: ±5% separation location

3. **Resolved turbulence content**:
   - Verify LES region has fluctuations (not RANS steady)
   - Check k_resolved / k_total > 0.8 in LES mode

**Acceptance Criteria**:
- ✅ Shielding function transitions correctly
- ✅ Separation prediction within ±10% of benchmark
- ✅ LES region shows resolved turbulence content

---

### Sprint 1.99.0 - Pressure Correction Derivation Documentation (2h)

**Objective**: Provide rigorous mathematical foundation for SIMPLE/PISO algorithms

#### Task 2.6: Document Pressure Correction Equation Derivation (2h)

**Current Gap**: Poisson equation stated without derivation

**Complete Derivation** (to be added in `pressure.rs:1-70`):

```rust
//! # Pressure Correction Equation Derivation
//!
//! ## Starting Point: Momentum Equation
//! Discretize momentum equation at time step n+1:
//! ```text
//! (ρV/Δt)(u*-un) = -∇p* + ν∇²u* + convection + sources
//! ```
//! where u* is predicted velocity (not divergence-free), p* is guessed pressure.
//!
//! ## Pressure Correction Concept
//! Define corrections: p = p* + p', u = u* + u'
//!
//! Subtract predicted from final momentum equation:
//! ```text
//! (ρV/Δt)(u' - 0) = -∇p' + ν∇²u' + ...
//! ```
//!
//! **SIMPLE Assumption**: Neglect ν∇²u' and neighbor velocity corrections
//! ```text
//! u' ≈ -(Δt/ρV)∇p'
//! ```
//!
//! ## Divergence-Free Constraint
//! Final velocity must satisfy continuity:
//! ```text
//! ∇·u^{n+1} = ∇·(u* + u') = 0
//! ```
//!
//! Substitute u' approximation:
//! ```text
//! ∇·u* - ∇·[(Δt/ρV)∇p'] = 0
//! ∇·[(Δt/ρV)∇p'] = ∇·u*
//! ```
//!
//! ## Poisson Equation for Pressure Correction
//! Assuming constant Δt and ρ:
//! ```text
//! ∇²p' = (ρ/Δt)∇·u*
//! ```
//! This is the **pressure correction Poisson equation** (Eq. 6.3-1 Patankar 1980).
//!
//! ## Boundary Conditions
//! - **Solid walls**: ∂p'/∂n = 0 (Neumann BC, no normal velocity correction)
//! - **Inlets**: p' = 0 if velocity specified (Dirichlet)
//! - **Outlets**: ∂p'/∂n = 0 or p' = 0 depending on BC type
//!
//! ## Reference Pressure
//! Poisson equation with pure Neumann BCs yields solution up to constant.
//! Fix one point: p'(i,j) = 0 OR use integral constraint ∫p' dV = 0.
//!
//! ## References
//! - Patankar, S.V. (1980). *Numerical Heat Transfer and Fluid Flow*. §6.3
//! - Ferziger & Perić (2002). *Computational Methods for Fluid Dynamics*. §7.5
//! - Issa, R.I. (1986). "Solution of the implicitly discretised fluid flow equations by operator-splitting." J. Comput. Phys. 62(1):40-65.
```

**Acceptance Criteria**:
- ✅ Complete step-by-step derivation in documentation
- ✅ Literature citations with specific equation numbers
- ✅ Boundary condition justifications included
- ✅ Reference pressure handling explained

---

## Phase 3: Documentation & Theory Enhancement (Sprints 1.100.0-1.103.0, 12h)

### Sprint 1.100.0 - Linear Solver Convergence Theory (6h)

**Objective**: Document convergence theorems, bounds, and optimality proofs

#### Task 3.1: GMRES Convergence Theory Documentation (2h)

**Add to** `crates/cfd-math/src/linear_solver/gmres/mod.rs`:

```rust
//! ## Convergence Theory
//!
//! ### Minimal Residual Property
//! GMRES minimizes ||b - Axₖ|| over x₀ + Kₘ(A,r₀), where Kₘ is the m-dimensional
//! Krylov subspace:
//! ```text
//! Kₘ(A,r₀) = span{r₀, Ar₀, A²r₀, ..., A^{m-1}r₀}
//! ```
//!
//! **Theorem (Saad & Schultz 1986)**: For any x₀, GMRES(m) finds xₘ such that:
//! ```text
//! ||rₘ|| = min_{x ∈ x₀ + Kₘ} ||b - Ax||
//! ```
//!
//! ### Convergence Bound
//! **Field of Values Bound** (Elman 1982):
//! ```text
//! ||rₘ|| / ||r₀|| ≤ κ(V) · min_{pₘ ∈ Πₘ, pₘ(0)=1} max_{z ∈ W(A)} |pₘ(z)|
//! ```
//! where W(A) is the field of values (numerical range) of A, Πₘ is polynomials of degree m.
//!
//! ### Practical Implications
//! - **Clustered eigenvalues**: Superlinear convergence expected
//! - **Outlier eigenvalues**: May slow convergence; preconditioning recommended
//! - **Restart effects**: GMRES(m) may stagnate if m too small; typical m ∈ [20,50]
//!
//! ### Complexity
//! - **Per iteration**: O(nnz) + O(m²) for Arnoldi orthogonalization
//! - **Memory**: O(n·m) for Krylov basis + O(m²) for Hessenberg matrix
//! - **Breakdown**: If ||vₖ|| < ε (happy breakdown), exact solution found in Kₘ
//!
//! ## References
//! - Saad, Y., & Schultz, M. H. (1986). GMRES: A generalized minimal residual algorithm.
//!   SIAM J. Sci. Stat. Comput., 7(3), 856-869. [Eq. 2.12 minimal residual property]
//! - Elman, H. C. (1982). Iterative methods for large sparse nonsymmetric systems.
//!   Yale University Tech Report 229. [Field of values bound]
//! - Saad, Y. (2003). Iterative Methods for Sparse Linear Systems (2nd ed.). §6.5
```

**Acceptance Criteria**:
- ✅ Minimal residual theorem stated with proof reference
- ✅ Convergence bound documented (Field of Values)
- ✅ Restart effects explained
- ✅ Complexity analysis included

---

#### Task 3.2: CG Optimality Theorem Documentation (2h)

**Add to** `crates/cfd-math/src/linear_solver/conjugate_gradient.rs`:

```rust
//! ## Mathematical Foundation
//!
//! ### Optimality Theorem
//! **Theorem (Hestenes & Stiefel 1952)**: Conjugate Gradient minimizes the A-norm
//! error ||xₖ - x*||_A over x₀ + Kₖ(A,r₀):
//! ```text
//! ||xₖ - x*||_A = min_{x ∈ x₀ + Kₖ} ||x - x*||_A
//! ```
//! where ||e||_A = √(eᵀAe) is the energy norm for SPD matrix A.
//!
//! ### Finite Termination
//! **Corollary**: In exact arithmetic, CG converges in ≤n iterations for n×n matrix A.
//!
//! **Proof Sketch**: Krylov subspace Kₙ(A,r₀) = ℝⁿ for generic r₀, so minimization
//! over Kₙ recovers exact solution x*.
//!
//! ### Convergence Rate
//! **Theorem (Chebyshev Polynomial Bound)**:
//! ```text
//! ||eₖ||_A / ||e₀||_A ≤ 2 [(√κ - 1) / (√κ + 1)]^k
//! ```
//! where κ = κ(A) = λₘₐₓ/λₘᵢₙ is the condition number of A.
//!
//! **Implication**: CG converges in O(√κ) iterations to reduce error by constant factor.
//!
//! ### Preconditioning Effect
//! Preconditioned CG solves M⁻¹Ax = M⁻¹b with convergence rate:
//! ```text
//! O(√κ(M⁻¹A))
//! ```
//! **Goal**: Choose preconditioner M ≈ A such that κ(M⁻¹A) << κ(A).
//!
//! ### Breakdown Prevention
//! CG requires A symmetric positive definite (SPD):
//! - **Symmetry**: Aᵀ = A ensures conjugacy relation
//! - **Positive definite**: xᵀAx > 0 for all x≠0 guarantees α = rᵀz / pᵀAp > 0
//!
//! If pᵀAp ≤ 0, matrix is indefinite → use BiCGSTAB or GMRES instead.
//!
//! ## References
//! - Hestenes, M. R., & Stiefel, E. (1952). Methods of conjugate gradients.
//!   J. Res. Nat. Bur. Stand., 49(6), 409-436. [Original CG algorithm]
//! - Shewchuk, J. R. (1994). An Introduction to the Conjugate Gradient Method.
//!   Carnegie Mellon Tech Report. [Excellent tutorial with proofs]
//! - Saad, Y. (2003). Iterative Methods for Sparse Linear Systems. §6.7
```

**Acceptance Criteria**:
- ✅ Optimality theorem stated with A-norm minimization
- ✅ Finite termination proven
- ✅ Convergence rate O(√κ) derived
- ✅ Preconditioning effect explained

---

#### Task 3.3: BiCGSTAB Literature References (2h)

**Add to** `crates/cfd-math/src/linear_solver/bicgstab.rs`:

```rust
//! ## Algorithm Overview
//!
//! BiCGSTAB (Biconjugate Gradient Stabilized) combines the conjugacy property of BiCG
//! with GMRES-like minimization to improve convergence for nonsymmetric systems.
//!
//! ### Why BiCGSTAB vs BiCG?
//! **BiCG Issue**: Irregular convergence with large oscillations in residual norm.
//!
//! **BiCGSTAB Solution**: Applies local GMRES(1) minimization at each iteration:
//! ```text
//! rₖ₊₁ = s - ω·t, where ω = argmin||s - ω·t||
//! ```
//! This stabilizes convergence while maintaining low memory O(n).
//!
//! ### Convergence Characteristics
//! - **Smoothness**: Residual norm decreases more steadily than BiCG
//! - **Speed**: Often faster than BiCG, competitive with GMRES(m) for moderate m
//! - **Memory**: Constant 8 vectors vs GMRES(m) requiring m+1 vectors
//!
//! ### Breakdown Conditions
//! 1. **Primary breakdown**: ρₖ = r̂ᵀrₖ ≈ 0 (shadow residual orthogonality loss)
//!    - Occurs when r̂ ⊥ Kₖ(A,r₀)
//!    - Remedy: Restart with different r̂ OR switch to GMRES
//!
//! 2. **Secondary breakdown**: ωₖ ≈ 0 (stabilization failure)
//!    - Occurs when s ≈ 0 but t ≠ 0
//!    - Remedy: Set ωₖ = 1 (reverts to BiCG step)
//!
//! ### When to Use BiCGSTAB
//! - **Nonsymmetric matrices**: Convection-diffusion, Navier-Stokes
//! - **Memory constrained**: Large problems where GMRES(m) storage prohibitive
//! - **Moderate condition number**: κ(A) < 10⁴ (higher may stagnate)
//!
//! ## References
//! - Van der Vorst, H. A. (1992). "Bi-CGSTAB: A fast and smoothly converging variant
//!   of Bi-CG for the solution of nonsymmetric linear systems."
//!   SIAM J. Sci. Stat. Comput., 13(2), 631-644. [Original algorithm]
//! - Sleijpen, G. L. G., & Fokkema, D. R. (1993). "BiCGstab(ℓ) for linear equations
//!   involving unsymmetric matrices with complex spectrum."
//!   ETNA, 1, 11-32. [BiCGSTAB(ℓ) variants]
//! - Saad, Y. (2003). Iterative Methods for Sparse Linear Systems. §7.4
```

**Acceptance Criteria**:
- ✅ Van der Vorst (1992) cited as primary reference
- ✅ BiCG vs BiCGSTAB differences explained
- ✅ Breakdown conditions documented with remedies
- ✅ Usage guidelines provided

---

### Sprint 1.101.0 - AMG Convergence Theory (4h)

**Objective**: Document algebraic multigrid theoretical foundations

#### Task 3.4: AMG Two-Grid Convergence Analysis (2h)

**Add to** `crates/cfd-math/src/linear_solver/multigrid/mod.rs`:

```rust
//! ## Algebraic Multigrid (AMG) Convergence Theory
//!
//! ### Two-Grid Method
//! A single V-cycle iteration consists of:
//! 1. **Pre-smoothing**: ν₁ iterations of smoother (e.g., Gauss-Seidel)
//! 2. **Restriction**: Compute coarse residual rᶜ = R(rᶠ)
//! 3. **Coarse solve**: Solve Aᶜeᶜ = rᶜ
//! 4. **Prolongation**: Interpolate correction eᶠ = P(eᶜ)
//! 5. **Post-smoothing**: ν₂ smoother iterations
//!
//! ### Convergence Theorem
//! **Theorem (Ruge & Stüben 1987)**: For SPD matrices, the two-grid convergence
//! factor ρ₂ satisfies:
//! ```text
//! ρ₂ ≤ max{||Sⁿᵘ¹||_A, ||Sⁿᵘ²(I - PAᶜ⁻¹RAᶠ)||_A}
//! ```
//! where S is the smoothing iteration matrix, ν₁,ν₂ are smoothing counts,
//! P is prolongation, R is restriction.
//!
//! **Optimal AMG**: Achieve ρ₂ < 0.2 (factor 5 error reduction per V-cycle).
//!
//! ### Key Requirements
//! 1. **Smoothing Property**: S damps high-frequency error
//!    - Gauss-Seidel, Jacobi typically satisfy ||Sᵏeₕ||_A → 0 for oscillatory eₕ
//!
//! 2. **Approximation Property**: Coarse grid captures smooth error
//!    - Interpolation P should satisfy ||eₛ - Peᶜ||_A ≈ 0 for smooth eₛ
//!    - Ruge-Stüben coarsening selects points ensuring this property
//!
//! ### Strength of Connection
//! **Definition**: i strongly connects to j if:
//! ```text
//! |aᵢⱼ| ≥ θ · max_{k≠i} |aᵢₖ|
//! ```
//! where θ ∈ [0,1] is strength threshold (typically 0.25).
//!
//! **Interpretation**: Strong connections indicate error components that must be
//! represented on coarse grid.
//!
//! ### Interpolation Formula
//! **Classical AMG Interpolation** (Ruge-Stüben):
//! ```text
//! (Peᶜ)ᵢ = {
//!     eᶜᵢ                           if i ∈ C (coarse point)
//!     -Σⱼ∈Cᵢ (aᵢⱼ/aᵢᵢ) eᶜⱼ          if i ∈ F (fine point)
//! }
//! ```
//! where Cᵢ is set of coarse points strongly influencing i.
//!
//! ### Multigrid Convergence
//! **Theorem**: For W-cycle with sufficient coarse levels, convergence factor:
//! ```text
//! ρ_MG ≤ C · ρ₂^L
//! ```
//! where L is number of levels, C is a small constant.
//!
//! **Implication**: O(L) = O(log n) work per iteration with exponentially small
//! error reduction ρ_MG ≈ 0.1-0.2.
//!
//! ## References
//! - Ruge, J. W., & Stüben, K. (1987). "Algebraic multigrid."
//!   In Multigrid Methods (pp. 73-130). SIAM. [Classical AMG theory]
//! - Briggs, W. L., Henson, V. E., & McCormick, S. F. (2000).
//!   A Multigrid Tutorial (2nd ed.). SIAM. [Accessible introduction]
//! - Trottenberg, U., Oosterlee, C. W., & Schüller, A. (2001).
//!   Multigrid. Academic Press. [Comprehensive theory]
//! - Stüben, K. (2001). "A review of algebraic multigrid."
//!   J. Comput. Appl. Math., 128(1-2), 281-309. [Modern AMG variants]
```

**Acceptance Criteria**:
- ✅ Two-grid convergence theorem stated with bound
- ✅ Smoothing and approximation properties explained
- ✅ Strength of connection defined
- ✅ Interpolation formula documented
- ✅ Multigrid convergence rate analyzed

---

#### Task 3.5: AMG Parameter Sensitivity Documentation (2h)

**Add parameter tuning guide** to `multigrid/mod.rs`:

```rust
//! ## Parameter Tuning Guide
//!
//! ### Strength Threshold (θ)
//! **Default**: 0.25  
//! **Range**: [0.0, 1.0]
//!
//! **Effect**:
//! - **θ → 0**: More strong connections → larger coarse grid → higher setup cost
//! - **θ → 1**: Fewer strong connections → smaller coarse grid → weaker interpolation
//!
//! **Tuning**:
//! - **Anisotropic problems** (e.g., boundary layers): θ ∈ [0.0, 0.25]
//! - **Isotropic problems** (e.g., Poisson): θ ∈ [0.25, 0.5]
//! - Monitor operator complexity: target Cₒₚ ∈ [1.2, 1.8]
//!
//! ### Smoothing Iterations (ν₁, ν₂)
//! **Default**: ν₁=2 (pre), ν₂=2 (post)  
//! **Range**: [1, 5]
//!
//! **Effect**:
//! - More smoothing → better high-frequency damping → slower per cycle
//! - Typical: ν₁ + ν₂ = 3 or 4
//!
//! **Tuning**:
//! - **Symmetric problems**: ν₁ = ν₂ (balanced)
//! - **Nonsymmetric problems**: ν₂ > ν₁ (more post-smoothing)
//!
//! ### Coarsening Strategy
//! **Options**: Ruge-Stüben, Aggregation, Hybrid
//!
//! **Ruge-Stüben**:
//! - Best for structured grids, anisotropic problems
//! - Higher setup cost O(nnz)
//! - Target: 2-4x coarsening ratio
//!
//! **Aggregation**:
//! - Faster setup, suitable for unstructured meshes
//! - May require more V-cycles to converge
//! - Target: 8-27x coarsening ratio (block-based)
//!
//! **Hybrid** (recommended):
//! - Ruge-Stüben on fine levels
//! - Aggregation on coarse levels
//! - Balances setup cost and convergence
//!
//! ### Cycle Type
//! **V-cycle**: γ=1 (visits each level once)  
//! **W-cycle**: γ=2 (visits each level twice)  
//! **F-cycle**: Hybrid between V and W
//!
//! **Trade-off**:
//! - V-cycle: O(n) per iteration, ρ ≈ 0.2-0.3
//! - W-cycle: O(n log n) per iteration, ρ ≈ 0.1-0.15
//!
//! **Recommendation**: Start with V-cycle; use W-cycle if convergence slow.
```

**Acceptance Criteria**:
- ✅ Strength threshold tuning guidelines
- ✅ Smoothing iteration recommendations
- ✅ Coarsening strategy comparisons
- ✅ Cycle type trade-offs documented

---

### Sprint 1.102.0 - Under-Relaxation Theory (2h)

**Objective**: Provide stability analysis for SIMPLE/PISO under-relaxation parameters

#### Task 3.6: Under-Relaxation Stability Derivation (2h)

**Add to** `crates/cfd-2d/src/simplec_pimple/config.rs`:

```rust
//! ## Under-Relaxation Theory
//!
//! ### Purpose
//! Under-relaxation stabilizes iterative algorithms for nonlinear coupled equations
//! by limiting solution changes per iteration:
//! ```text
//! φⁿ⁺¹ = φⁿ + α(φ* - φⁿ)
//! ```
//! where φ* is updated value from discrete equations, α ∈ (0,1] is relaxation factor.
//!
//! ### SIMPLE Algorithm Stability
//! **Theorem (Patankar 1980, §6.7)**: For SIMPLE algorithm, stability requires:
//! ```text
//! αᵤ ∈ (0, 1],  αₚ ∈ (0, 1 - αᵤ/2]
//! ```
//! where αᵤ is velocity relaxation, αₚ is pressure relaxation.
//!
//! **Derivation Sketch**:
//! 1. Momentum equation iteration matrix eigenvalues bounded by (1-αᵤ)
//! 2. Pressure correction coupling introduces αₚ·αᵤ cross-term
//! 3. Stability: spectral radius ρ(M) < 1 requires αₚ < 1 - αᵤ/2
//!
//! ### Practical Guidelines
//! **Standard Choices**:
//! - Laminar flows (Re < 100): αᵤ=0.7, αₚ=0.3
//! - Transitional (100 < Re < 1000): αᵤ=0.5-0.7, αₚ=0.2-0.3
//! - Turbulent (Re > 1000): αᵤ=0.4-0.6, αₚ=0.1-0.2
//!
//! **Why αₚ < αᵤ?**
//! Pressure correction depends on velocity field; slower pressure updates allow
//! velocity to stabilize first, reducing oscillations.
//!
//! ### SIMPLEC Modification
//! **SIMPLEC** (Van Doormaal & Raithby 1984) allows larger αₚ:
//! ```text
//! αᵤ ∈ (0,1],  αₚ ∈ (0, 1]  (less restrictive)
//! ```
//!
//! **Mechanism**: SIMPLEC modifies momentum coefficient aₚ → aₚ - Σaₙᵦ, reducing
//! neighbor influence and improving pressure-velocity coupling.
//!
//! **Typical SIMPLEC**: αᵤ=0.7, αₚ=0.5-0.7 (higher than SIMPLE)
//!
//! ### PIMPLE Relaxation
//! **PIMPLE** uses transient term for stability, allowing αₚ ≈ 1 if Δt small:
//! ```text
//! Stability: CFL · αᵤ < 1
//! ```
//! where CFL = max|u|Δt/Δx.
//!
//! **Adaptive Strategy**:
//! - Outer iterations < 3: Use under-relaxation (αᵤ=0.7, αₚ=0.3)
//! - Outer iterations ≥ 3: Reduce relaxation (αᵤ=0.9, αₚ=0.5)
//!
//! ### Dynamic Adjustment
//! **Convergence-based tuning**:
//! ```rust
//! if residual_decrease < 0.1 {
//!     // Convergence stalled → reduce relaxation
//!     alpha_u *= 0.9;
//!     alpha_p *= 0.9;
//! } else if residual_decrease > 0.5 && iterations < 5 {
//!     // Fast convergence → increase relaxation
//!     alpha_u = min(1.0, alpha_u * 1.1);
//!     alpha_p = min(1.0, alpha_p * 1.1);
//! }
//! ```
//!
//! ## References
//! - Patankar, S. V. (1980). Numerical Heat Transfer and Fluid Flow. §6.7
//! - Van Doormaal, J. P., & Raithby, G. D. (1984). "Enhancements of the SIMPLE
//!   method for predicting incompressible fluid flows." Numer. Heat Transfer, 7(2), 147-163.
//! - Ferziger, J. H., & Perić, M. (2002). Computational Methods for Fluid Dynamics. §7.5
```

**Acceptance Criteria**:
- ✅ Patankar stability criterion derived
- ✅ αₚ < αᵤ relationship explained
- ✅ SIMPLEC vs SIMPLE differences documented
- ✅ PIMPLE CFL-based stability included
- ✅ Adaptive tuning strategy provided

---

## Phase 4: Extended Test Coverage (Sprints 1.104.0-1.106.0, 8h)

### Sprint 1.104.0 - Backward-Facing Step Benchmark (4h)

**Objective**: Implement standard turbulent flow separation benchmark

#### Task 4.1: Backward-Facing Step Implementation (3h)

**Benchmark Details**:
- **Geometry**: 2D channel with sudden expansion (height ratio 2:1)
- **Reynolds**: Re_h = 37,500 (based on step height h)
- **Reference**: Driver & Seegmiller (1985), Armaly et al. (1983)

**Setup**:
1. **Domain** (1h):
   - Upstream: -4h to 0
   - Step: height h
   - Downstream: 0 to 30h
   - Upper wall: 2h
   - Mesh: ~20,000 cells, y+=1 at walls

2. **Boundary Conditions** (30min):
   - Inlet: Fully developed turbulent profile (1/7 power law)
   - Outlet: Zero gradient (convective BC)
   - Walls: No-slip + wall functions

3. **Solver Settings** (30min):
   - Algorithm: SIMPLEC or PIMPLE
   - Turbulence: k-ω SST (recommended for separation)
   - Under-relaxation: αᵤ=0.7, αₚ=0.3
   - Convergence: Continuity residual < 1e-4

**Validation Metrics**:
- **Reattachment length** x_r/h (primary metric):
  - Driver & Seegmiller (1985): x_r/h = 6.0 ± 0.5
  - Armaly et al. (1983): x_r/h = 6.28
- **Velocity profiles**: u(y) at x/h = 2, 4, 6, 8
- **Skin friction**: Cƒ vs x/h along bottom wall

**Acceptance Criteria**:
- ✅ Reattachment length within ±10% of experimental: x_r/h ∈ [5.4, 6.6]
- ✅ Velocity profile L2 error < 15% vs Driver & Seegmiller data
- ✅ Skin friction shows negative region (separation bubble)

---

#### Task 4.2: Grid Independence Study (1h)

**Objective**: Verify solution independence from mesh resolution

**Grids**:
1. Coarse: ~5,000 cells
2. Medium: ~20,000 cells
3. Fine: ~80,000 cells

**Procedure**:
1. Run all three grids to convergence
2. Compare x_r/h for each mesh
3. Compute Grid Convergence Index (GCI) per Roache (1998)
4. Verify GCI < 5% (asymptotic range)

**Deliverable**:
- Table: Mesh size vs x_r/h
- Plot: x_r convergence with Richardson extrapolation
- GCI calculation confirming mesh independence

**Acceptance Criteria**:
- ✅ GCI < 5% indicates asymptotic range reached
- ✅ Medium and fine mesh results agree within 2%

---

### Sprint 1.105.0 - Rectangular Channel Extended Tests (2h)

**Objective**: Comprehensive validation for non-square geometries

#### Task 4.3: Aspect Ratio Dependency Tests (1h)

**Test Matrix**:
- Aspect ratios: α ∈ {0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0}
- Reynolds numbers: Re ∈ {50, 100, 500, 1000, 2000}
- Total: 35 test cases

**Validation**:
1. **Poiseuille Number Correlation** (Shah & London 1978):
   - Wide α ≥ 1: 5-term polynomial (lines 150-167)
   - Tall α ≤ 1: Alternative polynomial (lines 168-188)
   - Square α=1: Po = 56.91 exactly

2. **Velocity Profile**:
   - Compare with analytical solution u(y,z) from separation of variables
   - Check maximum velocity location shifts with aspect ratio

3. **Pressure Drop Scaling**:
   - Verify ΔP ∝ μ·L·V (linear dependencies)
   - Check ΔP ∝ 1/Dh² (geometric scaling)

**Acceptance Criteria**:
- ✅ Po values match Shah & London within ±2%
- ✅ All 35 test cases pass
- ✅ Velocity profiles agree with analytical solution (L2 error <5%)

---

#### Task 4.4: Reynolds Number Variation Tests (1h)

**After fixing hardcoded Re=100 bug**, validate over full range:

**Test Cases**:
- Low Re: {10, 50, 100}
- Moderate Re: {500, 1000, 1500}
- Transition: {2000, 2100, 2200, 2300}

**Validation**:
1. **Laminar Regime** (Re < 2300):
   - f·Re = Po (constant)
   - Verify this relationship holds

2. **Transition Detection**:
   - At Re=2300, should warn or error (exceeds validity)
   - Ensure model doesn't silently produce wrong results for Re > 2300

3. **Scaling Consistency**:
   - Doubling Re → halves f_fanning
   - Pressure drop ΔP ∝ Re for fixed geometry

**Acceptance Criteria**:
- ✅ f·Re = Po ± 1% for all Re < 2300
- ✅ Re > 2300 triggers validity warning
- ✅ Scaling laws verified

---

### Sprint 1.106.0 - AMG Parameter Sensitivity Study (2h)

**Objective**: Empirically validate AMG convergence theory predictions

#### Task 4.5: Strength Threshold Sensitivity (1h)

**Test Setup**:
- Problem: 2D Poisson equation on 256×256 grid (65,536 unknowns)
- Right-hand side: f(x,y) = sin(πx)sin(πy)
- Boundary: Dirichlet u=0

**Parameter Sweep**:
- θ ∈ {0.0, 0.1, 0.2, 0.25, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8}

**Metrics**:
1. **Setup time** vs θ
2. **Operator complexity** Cₒₚ vs θ
3. **Convergence factor** ρ vs θ
4. **Time to solution** (setup + iterations) vs θ

**Expected Behavior** (from theory):
- θ → 0: High Cₒₚ, low ρ (good convergence, expensive setup)
- θ → 1: Low Cₒₚ, high ρ (cheap setup, poor convergence)
- Optimal: θ ∈ [0.2, 0.3] for Poisson problems

**Deliverable**:
- 4 plots: Setup time, Cₒₚ, ρ, total time vs θ
- Table summarizing optimal θ ranges for different problem types

**Acceptance Criteria**:
- ✅ Optimal θ ≈ 0.25 ± 0.05 confirms theory
- ✅ Operator complexity Cₒₚ ∈ [1.2, 1.8] at θ=0.25
- ✅ Convergence factor ρ < 0.3 at θ=0.25

---

#### Task 4.6: Cycle Type Comparison (1h)

**Comparison**: V-cycle vs W-cycle performance

**Test Problems**:
1. Isotropic Poisson (64² grid)
2. Anisotropic diffusion (aspect ratio 100:1, 64² grid)
3. Convection-diffusion (Pe=100, 64² grid)

**Metrics**:
- Iterations to convergence (residual < 1e-6)
- Work per cycle (SPMV count)
- Total work = iterations × work/cycle

**Theory Predictions**:
- V-cycle: ρ ≈ 0.2-0.3, O(n) work per cycle
- W-cycle: ρ ≈ 0.1-0.15, O(n log n) work per cycle

**Deliverable**:
- Table: V vs W for 3 test problems
- Analysis: When is W-cycle worth the extra cost?

**Acceptance Criteria**:
- ✅ V-cycle 20-30% faster for isotropic problems
- ✅ W-cycle 10-20% faster for anisotropic/difficult problems
- ✅ Results match theoretical complexity predictions

---

## Sprint-by-Sprint Summary

| Sprint | Phase | Focus Area | Hours | Critical Issues Resolved |
|--------|-------|------------|-------|-------------------------|
| 1.95.0 | 1 | Critical Bugs | 4h | Rectangular Re=100 hardcoding |
| 1.96.0 | 1 | SIMPLE Accuracy | 4h | Root cause analysis for 30% errors |
| 1.97.0 | 2 | Entrance Effects | 6h | Theory, tests, coefficient justification |
| 1.98.0 | 2 | DES Enhancement | 4h | Proper Spalart 2006 shielding function |
| 1.99.0 | 2 | Pressure Derivation | 2h | Complete Poisson equation derivation |
| 1.100.0 | 3 | Linear Solver Theory | 6h | GMRES, CG, BiCGSTAB convergence proofs |
| 1.101.0 | 3 | AMG Theory | 4h | Two-grid convergence, parameter tuning |
| 1.102.0 | 3 | Under-Relaxation | 2h | SIMPLE/SIMPLEC stability analysis |
| 1.104.0 | 4 | Backward-Facing Step | 4h | Standard separation benchmark |
| 1.105.0 | 4 | Rectangular Tests | 2h | Aspect ratio, Reynolds validation |
| 1.106.0 | 4 | AMG Sensitivity | 2h | Empirical parameter validation |

**Total Effort**: 40 hours over 11 sprints

---

## Success Metrics

### Phase 1 Completion (Sprint 1.96.0)
- ✅ Zero critical bugs (P0 issues resolved)
- ✅ SIMPLE/SIMPLEC root cause identified
- ✅ Path forward for accuracy improvement defined

### Phase 2 Completion (Sprint 1.99.0)
- ✅ All high-priority gaps (P1) addressed
- ✅ Entrance effects: Theory + tests implemented
- ✅ DES: Production-grade shielding function
- ✅ SIMPLE: Complete mathematical derivation

### Phase 3 Completion (Sprint 1.102.0)
- ✅ All linear solvers: Convergence theory documented
- ✅ AMG: Two-grid analysis with parameter guidance
- ✅ SIMPLE: Under-relaxation stability criterion

### Phase 4 Completion (Sprint 1.106.0)
- ✅ Backward-facing step: Standard benchmark passing
- ✅ Rectangular channel: Comprehensive test coverage
- ✅ AMG: Empirical parameter validation
- ✅ **Target**: >85% mathematical rigor across all 17 audited components

---

## Risk Mitigation

### Risk 1: SIMPLE Accuracy Cannot Be Fixed (High Impact)

**Mitigation**:
- Sprint 1.96.0: Identify if issue is fundamental algorithm limitation
- If unfixable: Document limitation clearly in README
- Alternative: Recommend PISO or PIMPLE for production use
- Deprecation path: Mark SIMPLE as "experimental" if accuracy inadequate

### Risk 2: Entrance Effects Theory Unavailable (Medium Impact)

**Mitigation**:
- Sprint 1.97.0: If derivations not found, cite empirical sources
- Document uncertainty: "Correlation valid for Re ∈ [X,Y] per reference"
- Implement conservative error bounds
- Add validation against experimental data (Idelchik handbook)

### Risk 3: DES Benchmark Data Inaccessible (Low Impact)

**Mitigation**:
- Sprint 1.98.0: Use qualitative validation (shielding function profile)
- Check fd ∈ [0,1] in boundary layer vs freestream
- If no DNS data: Document limitation, proceed with implementation
- Future work: Validate when benchmark data becomes available

---

## Continuous Integration

### Automated Quality Gates (Enforce After Phase 2)

```yaml
# .github/workflows/mathematical_rigor.yml
name: Mathematical Rigor CI

on: [push, pull_request]

jobs:
  test-rigor:
    steps:
      - name: Check for hardcoded magic numbers
        run: |
          # Fail if new hardcoded values without documentation
          grep -n "let.*=.*[0-9]\+\.[0-9]\+" src/**/*.rs | \
          grep -v "// Justified:" && exit 1 || exit 0

      - name: Verify literature references
        run: |
          # Ensure all theorems cite sources
          grep -n "Theorem\|Equation" src/**/*.rs | \
          grep -v "Reference:\|Citation:" && exit 1 || exit 0

      - name: Test coverage check
        run: |
          # Minimum 80% coverage for mathematical components
          cargo tarpaulin --out Xml
          coverage=$(grep -oP 'line-rate="\K[0-9.]+' cobertura.xml)
          if (( $(echo "$coverage < 0.80" | bc -l) )); then
            echo "Coverage $coverage < 80%"; exit 1
          fi
```

---

## Conclusion

This roadmap systematically elevates CFDrs from **mixed rigor (50-100%)** to **production-grade completeness (>85%)** by:

1. **Immediate fixes** (Sprints 1.95-1.96): Resolving showstopper bugs
2. **Theory enhancement** (Sprints 1.97-1.99): Filling critical mathematical gaps
3. **Documentation** (Sprints 1.100-1.102): Providing convergence proofs and stability analysis
4. **Validation expansion** (Sprints 1.104-1.106): Implementing standard benchmarks

**Estimated Total Effort**: 40 hours over 11 sprints  
**Target Completion**: Sprint 1.106.0  
**Expected Outcome**: All 17 audited components achieve >85% mathematical rigor

With disciplined execution of this plan, CFDrs will reach **commercial-grade theoretical foundations** comparable to ANSYS Fluent and OpenFOAM.

---

**Document Status**: ROADMAP COMPLETE  
**Next Action**: Stakeholder review and sprint 1.95.0 initiation  
**Owner**: CFD Development Team  
**Priority**: HIGH (production blocker resolution)
