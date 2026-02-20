# Gap Audit: `cfd-1d` Mathematical Purity (CRITICAL-011)

## Overview
An audit of `cfd-1d` was conducted to ensure all algorithms are strictly implemented from first principles, heavily documented with mathematical theorems, and devoid of placeholders or approximations, in accordance with the project's elite architectural constraints.

## Findings: Mathematical Approximations

### 1. Womersley Flow Velocity Profile (`womersley.rs`)
- **Current State**: Uses asymptotic closed-form approximations for low ($\alpha < 1$) and high ($\alpha > 10$) Womersley numbers, with an arbitrary interpolated blend for the transitional regime.
- **Gap**: The true velocity profile requires evaluation of the complex Bessel function of the first kind, zero order: $J_0(z)$ where $z = i^{3/2} \alpha \xi$. The current implementation fundamentally fails the correctness principle by approximating this complex arithmetic. 
- **Required Fix**: Implement an exact Bessel function calculator using power series for complex arguments (or Kelvin functions `ber`, `bei`) and re-derive the Womersley profile directly from the exact analytical solution.

### 2. Darcy Friction Factor (`darcy_weisbach.rs` & `channel/solver.rs`)
- **Current State**: Uses the explicit Haaland equation to approximate the Darcy friction factor for turbulent flow.
- **Gap**: The Haaland equation is only an approximation (error up to $\pm 2\%$) of the implicit Colebrook-White equation.
- **Required Fix**: The exact Colebrook equation ($1/\sqrt{f} = -2 \log_{10}(\epsilon/3.7D_h + 2.51/(Re\sqrt{f}))$) must be solved exactly via an iterative numerical method (e.g., Newton-Raphson) with rigorous convergence bounds.

### 3. Ellipse Perimeter / Hydraulic Diameter (`geometry.rs` & `branching/physics.rs`)
- **Current State**: Implements a "Ramanujan-like" approximation ($D_h \approx \sqrt{ab}$ or $2ab/(a+b)$) for the perimeter and hydraulic diameter of an ellipse.
- **Gap**: The perimeter of an ellipse is not algebraically solvable; it requires the Complete Elliptic Integral of the Second Kind, $E(e)$.
- **Required Fix**: Implement $E(e)$ exactly using the infinite series or AGM (Arithmetic-Geometric Mean) method to guarantee mathematically perfect boundary evaluation.

### 4. Serpentine Flow Approximations (`resistance/models/serpentine.rs`)
- **Current State**: Uses "White 1929 approximation" for negligible curvature effect and a "circular approximation" for shear rate estimation.
- **Gap**: Dean flow physics must be accurately accounted for using the exact Dean Number ($De$) rather than rough historical geometry simplifications.
- **Required Fix**: Apply exact Dean vortices physics models to precisely capture secondary flow resistances in rigid bends.

## Gap Audit: `cfd-mesh` Robustness and Mathematical Purity (CRITICAL-012)

### Overview
An audit of `cfd-mesh` was conducted focusing on Constructive Solid Geometry (CSG), BSP tree generation, vertex welding, and watertight snapping. The goal is to ensure all geometric processing aligns with the highest standards of mathematical correctness and eliminates heuristic-based approximations, based on the latest 2024 research in robust mesh booleans (e.g., EMBER: Exact Mesh Booleans via plane-based representations).

### Findings: Geometric Heuristics and Precision Gaps

#### 1. CSG Boolean Operations (`csg/boolean.rs`)
- **Current State**: Implements standard mesh face-soup booleans reliant on a heuristic containment check and naive BSP clipping. Uses `is_curved_mesh` to classify curved vs. flat surfaces by quantizing normals to a coarse grid (0.1 resolution).
- **Gap**: The containment detection and splitting logic suffer from floating-point inaccuracies, leading to sliver triangles, non-manifold edges, and shattering. The normal quantization is an unacceptable heuristic approximation. State-of-the-art methods (e.g., EMBER) mandate exact predicate evaluations or plane-based integer representations to guarantee exact, reliable booleans without degenerate failures.
- **Required Fix**: Refactor CSG booleans to use robust multi-precision exact predicates (like Shewchuk's or arithmetic expansions) and a plane-based representation. Eliminate normal quantization heuristics entirely.

#### 2. BSP Tree Splitting (`csg/bsp.rs`)
- **Current State**: Calculates the splitting plane using a greedy heuristic (`4 * spans + |front - back|`) over up to 16 candidates. Relies on `BSP_PLANE_TOLERANCE = 1e-5` for classification.
- **Gap**: The floating-point tolerance approach fails for coplanar or nearly coplanar intersections, violating mathematical purity. A heuristic candidate subset does not guarantee optimal or failure-free trees.
- **Required Fix**: Implement exact geometry predicates for BSP classification. Transition from point-based floating-point geometry to plane-based geometry (representing vertices as intersections of planes) to eliminate epsilon tolerances during inside/outside classification.

#### 3. Vertex Welding & Snapping (`welding/snap.rs` & `welding/welder.rs`)
- **Current State**: Uses an i64-quantized `SnappingGrid` with a 27-neighborhood search, bounded by a user-defined floating-point tolerance `eps` (e.g., 1e-6). `welder.rs` uses a basic `SpatialHashGrid`.
- **Gap**: While spatial hashing is fast, epsilon-welding without topological awareness breaks manifoldness (e.g., collapsing distinct but close surface sheets). True epsilon-robustness must prevent topological inversions.
- **Required Fix**: Upgrade from naive distance-based collapsing to topology-preserving vertex welding. Implement robust vertex merging that guarantees preservation of 2-manifold properties and avoids non-manifold edge creation, utilizing topological invariant checks during the merge phase.

#### 4. Curved Mesh Arrangement Pipeline (`csg/arrangement.rs`)
- **Current State**: Uses a 5-phase Boolean layout with exact algebraic coplanarity grouping via rigorous `orient3d`, strict GWN raycasting algorithm for inclusion, and mathematically robust barycentric evaluation.
- **Gap**: Fixed.
- **Status**: ✅ **FIXED** - Eliminated all $10^{-N}$ epsilon fallbacks and normal scaling approximations. Replaced heuristic ray-casting with exact mathematical algorithms.
