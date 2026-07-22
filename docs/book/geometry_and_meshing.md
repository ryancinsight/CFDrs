# Chapter 6 — Geometry Construction and CFD Coupling

## Overview

CFDrs geometry is built **constructively**: a **CSG primitive catalogue**
composes via Boolean operations with a **signed distance field (SDF)**
representation, feeds a **mesh generator**, and bonds to a **fluid-solid
mapping**. The result is a single `Mesh` or structured grid that any 1-D, 2-D,
or 3-D solver consumes. This chapter covers the primitive catalogue, Boolean
CSG algebra with SDF semantics, Delaunay/tetrahedral mesh generation, mesh
quality metrics and untangling, fluid-solid boundary wiring, the dimension-lifting
ladder (1-D network → 2-D channel → 3-D FEM), and the Atlas geometry types that
replace `nalgebra`.

---

## 1 — CSG Primitives

### Catalogue

`cfd-schematics::primitives` exposes the building blocks:

```rust
pub fn box_<F>(origin: Point3<F>, extents: Vector3<F>) -> Mesh;
pub fn cylinder<F>(origin: Point3<F>, axis: Vector3<F>, r: F, h: F) -> Mesh;
pub fn cone<F>(origin: Point3<F>, axis: Vector3<F>, r_base: F, r_tip: F, h: F) -> Mesh;
pub fn sphere_section<F>(origin: Point3<F>, r: F, theta_min: F, theta_max: F) -> Mesh;
pub fn torus<F>(center: Point3<F>, r_major: F, r_minor: F) -> Mesh;
pub fn capsule<F>(a: Point3<F>, b: Point3<F>, r: F) -> Mesh;
```

Additional domain-specific primitives in `cfd-schematics::bifurcation`,
`::venturi`, `::serpentine`, `::tpms`:

| Builder | Underlying CSG | Diagram |
|---|---|---|
| `BifurcationBuilder` | Two cylinders + fillet torus + optional cone stenosis | Y-shape with stenosis tag |
| `VenturiBuilder` | Two frustum cones + cylindrical throat | Converge-diverge axisymmetric channel |
| `SerpentineBuilder` | Chain of capsules (straight + half-torus turns) | Serpentine / Dean-vortex chip |
| `TpmsScaffold::gyroid` | Implicit SDF $F(x,y,z)=0$ offset by $t$ (porosity control) | Schwarz P / Gyroid |

Each primitive returns a `Mesh` whose vertices are `leto::Point3<F>` (Atlas) or
`nalgebra::Point3<F>` (legacy) and whose elements are 2-D / 3-D references into
the same `(i, j, k)` index space as the solver's grid.

### Atlas Geometry Types

The Atlas migration replaces `nalgebra` point/vector types with `leto`:

| Legacy | Atlas | Operation |
|---|---|---|
| `nalgebra::Point3<F>` | `leto::Point3<F>` | `.coords`, `.distance(&other)` |
| `nalgebra::Vector3<F>` | `leto::Vector3<F>` | `dot`, `cross`, `normalize` |
| `nalgebra::Isometry3<F>` | `leto::Isometry3<F>` | `*` (composition), `inverse()` |

All CSG operations are posed as `Point3 → SDF` or `Mesh → Mesh` depending on
the backend flag:

```rust
pub trait Sdf<F> {
    fn sdf(&self, p: Point3<F>) -> F;           // signed distance
    fn normal(&self, p: Point3<F>) -> Vector3<F>; // ∇SDF / |∇SDF|
    fn contains(&self, p: Point3<F>) -> bool { self.sdf(p) <= F::zero() }
}
```

---

## 2 — SDF Representation

### Signed Distance Definition

For a closed solid $S \subset \mathbb{R}^3$, the signed distance field is

$$d(\mathbf{x}) = \begin{cases}
-\min_{\mathbf{y}\in\partial S} \|\mathbf{x}-\mathbf{y}\| & \mathbf{x}\in\mathrm{int}(S)\\
+\min_{\mathbf{y}\in\partial S} \|\mathbf{x}-\mathbf{y}\| & \mathbf{x}\notin S
\end{cases} \tag{1}$$

with $|\nabla d| = 1$ almost everywhere (eikonal equation). SDF provides:

- **Inside-outside test**: $\mathrm{sign}(d(\mathbf{x}))$ in $O(1)$ without ray
  casting.
- **Normal**: $\mathbf{n}(\mathbf{x}) = \nabla d(\mathbf{x})$ for $\mathbf{x}\in\partial S$.
- **Offsetting**: surface at distance $t$ is simply $d(\mathbf{x}) = t$.

CFDrs SDFs are implemented both as analytical expressions (sphere, box,
cylinder — `cfd-schematics::sdf::analytical`) and as discretized grids
(tetrahedral mesh SDF via `cfd-schematics::sdf::mesh_sdf` using a BVH for
acceleration).

### Primitive SDFs

Sphere radius $R$ at origin: $d(\mathbf{p}) = \|\mathbf{p}\| - R$.

Box half-extents $\mathbf{b}$: with $\mathbf{q} = |\mathbf{p}| - \mathbf{b}$,

$$d_{\mathrm{box}}(\mathbf{p}) = \min(\max(q_x,q_y,q_z), 0) + \| \max(\mathbf{q}, 0) \| \tag{2}$$

Cylinder (infinite axis $z$): $d_{\mathrm{cyl}} = \sqrt{x^2+y^2} - R$.

Capped cylinder height $h$: $d_{\mathrm{cap}} = \max(d_{\mathrm{cyl}}, |z|-h/2)$.

Gyroid TPMS (implicit surface $F=0$, not a true SDF but treated as SDF proxy):

$$F(x,y,z) = \sin x\cos y + \sin y\cos z + \sin z\cos x - t = 0 \tag{3}$$

with $t \in [-1,1]$ controlling fluid volume fraction. True distance is
approximated by $|F|/\|\nabla F\|$ (first-order correction), sufficient for
meshing as long as the correction error is tracked.

---

## 3 — CSG Boolean Operations

### Set Algebra

Given two solids $A, B$ as SDFs $d_A, d_B$, the Boolean results are

$$A \cap B: \quad d_{A \cap B}(\mathbf{x}) = \max(d_A(\mathbf{x}), d_B(\mathbf{x})) \tag{4}$$

$$A \cup B: \quad d_{A \cup B}(\mathbf{x}) = \min(d_A(\mathbf{x}), d_B(\mathbf{x})) \tag{5}$$

$$A \setminus B: \quad d_{A \setminus B}(\mathbf{x}) = \max(d_A(\mathbf{x}), \ -d_B(\mathbf{x})) \tag{6}$$

Negation: $\neg A: d_{\neg A} = -d_A$.

Equations (4)-(6) are **exact SDFs** of the Boolean result only at points whose
closest boundary feature comes from a single operand; at crease/curved-crease
regions the true SDF differs but the sign (inside/outside test) remains exact.
For high-quality normals at creases, CFDrs uses the "exact closest feature"
SDF when available, falling back to smooth union with a fillet radius $k$:

$$d_{\mathrm{smooth-union}} = -\frac{1}{k}\ln\left(e^{-k d_A} + e^{-k d_B}\right) \approx \min(d_A, d_B) - \frac{\ln 2}{k} \text{ at } d_A=d_B \tag{7}$$

where $k \to \infty$ recovers the hard min, and finite $k$ gives a radius-$O(1/k)$
fillet blending creases — useful to avoid sharp re-entrant corners that generate
singular shear.

### Mesh-Level CSG

On mesh operands ($\text{Triangle mesh } A$, mesh $B$), three backends exist:

| Backend | Feature flag | Method | Topology |
|---|---|---|---|
| `csgrs` | `csgrs_backend` | BSP-tree mesh bool | Exact, slow |
| `brepkit` | `brepkit_backend` | B-rep NURBS boolean | Curved boundaries |
| `hephaestus` | `hephaestus_backend` | GPU SDF ray-march → mesh extraction via Marching Cubes / Dual Contouring | SDF, opt-in |

The trait:

```rust
pub trait CsgOp {
    fn union(self, other: Mesh) -> Mesh;
    fn intersection(self, other: Mesh) -> Mesh;
    fn difference(self, other: Mesh) -> Mesh;
}
```

Operations are deferred; a CFG flag selects the backend kernel. For SDF operands,
union/intersection are $O(1)$ (min/max); extraction to a mesh happens only when
`build_mesh()` is called.

### Complexity note

Mesh-mesh Boolean with $n$ triangles per mesh can produce up to $O(n^2)$
triangles in the result (naïve pairing). CFDrs BSP kernel is
$O(n \log n + k)$ with $k$ the number of output intersections; GPU SDF path
is $O(m)$ for a uniform grid of $m$ voxels independent of input mesh count.

---

## 4 — Mesh Generation

### Structured Grids

```rust
let grid = StructuredGrid2D::<f64>::unit_square(41, 41)?;
let grid = StructuredGrid2D::<f64>::clustered_wall(
    41, 41, Clustering::Tanh{ beta: 1.5 }, Clustering::Uniform,
)?;
```

Clustering functions map uniform computational coordinates $\xi \in [0,1]$ to
physical $x(\xi)$:

- **`Tanh`**: $x(\xi) = 1 - \tanh(\beta(1-\xi)) / \tanh(\beta)$ clusters near
  $x=0$ (wall) with strength $\beta$.
- **`Geometric`**: successive ratio $r$, total extent $L$, $x_j = L(r^j-1)/(r^n-1)$.
- **`Cosine` / Chebyshev**: $x_j = \tfrac12(1-\cos(\pi j/(n-1)))$ clusters at both
  ends.
- **`Uniform`**: $x(\xi) = \xi$.

For wall-resolved LES/DNS channel flow, `Tanh` with $\beta\approx 2$ places
the first cell at $y^+_1 < 1$ with $\sim 64$ cells.

### Tetrahedral Unstructured Mesh (3-D)

```rust
pub trait MeshGenerator<F: FloatElement> {
    fn generate(&self, geometry: &Mesh, density: &DensitySpec) -> Result<Mesh, CfdError>;
}
```

The default generator produces **tetrahedral** P1 elements for 3-D geometries and
**triangular** P1 elements for 2-D. Implementation uses Delaunay refinement
(Chew's second algorithm) with quality guard. Density control via
`DensitySpec`:

```rust
pub struct DensitySpec<F> {
    pub base_size: F,               // uniform background size
    pub refinements: Vec<(Aabb<F>, F)>, // per-box size override
    pub curvature_refinement: Option<F>, // angular deviation threshold (rad)
}
```

### SDF → Mesh Extraction

For SDF-based primitives (TPMS, implicit surfaces), marching cubes extracts
the isosurface $d(\mathbf{x}) = t$ at resolution `DensitySpec::base_size`.
CFDrs wraps `mc_sdf` (or GPU dual contouring under `hephaestus`) with shared
vertex deduplication:

```
For each grid cube edge (a,b):
  if d(a)*d(b) < 0:
    emit vertex at lerp(a,b) = a + d(a)/(d(a)-d(b)) * (b-a)
For each face with sign change: emit triangle per MC case table
Merge coincident vertices (hash map on quantized position)
```

Discretization error scales as $O(h^2)$ for MC without cubical feature recovery
and $O(h)$ for normals; DC recovers sharp features ($O(h^2)$).

---

## 5 — Mesh Quality Metrics

Built into `cfd-schematics::quality` and `cfd-3d::mesh::quality`:

### Triangle Quality

- **Aspect ratio**: $\mathrm{AR} = L_{\max}/L_{\min}$ or $R/\rho$ (circumradius / inradius),
  ideal $=1$ (equilateral).
- **Min angle**: $\theta_{\min} = \min(\theta_1, \theta_2, \theta_3)$;
  Delaunay triangulations maximize $\min_e \theta_{\min}^{(e)}$.
- **Area / signed area**: $\tfrac12|(\mathbf{b}-\mathbf{a})\times(\mathbf{c}-\mathbf{a})|$; zero or
  negative means degenerate / inverted.

### Tetrahedron Quality

- **Dihedral angle**: $\cos\theta_{ij} = -\mathbf{n}_i\cdot\mathbf{n}_j$,
  $\theta_{ij}$ between face normals; $\theta_{\min}>5^\circ$–$10^\circ$ guards
  against slivers.
- **Radius ratio** $\rho = 3 r_{\mathrm{in}}/r_{\mathrm{circ}}$ ∈ $[0,1]$,
  1 for regular tet; CFDrs rejects tets with $\rho < 0.1$.
- **Volume**: signed, $\mathrm{vol} = (\mathbf{b}-\mathbf{a})\cdot((\mathbf{c}-\mathbf{a})\times(\mathbf{d}-\mathbf{a}))/6$;
  nonpositive means inversion.

### Solver Relevance

| Metric | Impact when poor | Guard |
|---|---|---|
| Small dihedral angle (sliver) | Matrix conditioning $\kappa(A)\sim 1/\sin\theta_{\min}$ degrades Krylov convergence | Swap 4→4 or 5→6 flips, Laplacian smoothing |
| High aspect ratio | Diffusion discretization anisotropy, CFL penalty | Anisotropic refinement aware of flow direction |
| Inverted elements | Determinant sign flip → wrong advection sense | Delaunay re-triangulation repair |
| Non-orthogonal mesh (FV) | Pressure-momentum coupling convergence slowdown for SIMPLEC | `n_non_orthogonal_correctors` in PIMPLE loop |

### Smoothing and Optimization

- **Laplacian smoothing**: $\mathbf{x}_i^{\mathrm{new}} = (1-\omega)\mathbf{x}_i + \omega \frac{1}{|N(i)|}\sum_{j\in N(i)}\mathbf{x}_j$.
- **Optimization-based**: move vertices to maximize local quality measure (e.g. mean radius ratio) via L-BFGS per patch — used only for TPMS scaffold healing where Laplacian would collapse thin ligaments.

---

## 6 — Fluid-Solid Interface and Boundary Wiring

A fluid-solid interface is an `Interface` with two sides:

```rust
pub struct Interface<F: FloatElement> {
    fluid_id: RegionId,
    solid_id: RegionId,
    boundary_kind: BoundaryKind<F>,
}
```

The solver indexes its degree of freedom at the interface and applies the
chosen boundary condition:

- **No-slip** (`Wall`): $\mathbf{u}=0$ over solid surfaces.
- **Lid** (`Lid{ velocity }`): driven lid — $U_{\mathrm{lid}}$ tangent, $0$ normal.
- **Prescribed inlet** (`Inlet { parabolic { umax } }`): Hagen-Poiseuille
  parabolic used for pipe and bifurcation parent inlets.
- **Partial-slip** (`PartialSlip{ slip_length }`): relevant for TPMS scaffold
  surface roughness and hydrophobic channels.

For immersed boundary (IBM) geometries (`cfd-2d::ibm`), solid cells are masked
via a boolean array; BCs are applied at the fluid-solid cut including cut-cell
volume correction per Udaykumar et al. (2001).

---

## 7 — Dimension Lifting — 1-D to 3-D Fidelity Ladder

### Scenario Taxonomy

CFDrs examples operate at three fidelity levels linked by the same
`NetworkBlueprint` / `VenturiSpec` geometry description:

```
NetworkBlueprint  ──►  cfd-1d solve (network, fast)  ──►  AnalysisOverlay
        │
        ├──────────►  cfd-2d solve (cross-section, medium)
        │
        └──────────►  cfd-3d solve (full FEM/spectral, expensive)
        │
        └──────────►  fidelity comparison report
```

`cfd-1d` uses the hydraulic network abstraction:
$\Delta p_j = R_j Q_j + L_j \dot{Q}_j$ (resistance + inertance) with
Poiseuille-based $R$ and Womersley-based $L$. `cfd-2d` solves the full
Navier-Stokes in cross-section with immersed-boundary masking. `cfd-3d`
solves FEM Stokes or full spectral N-S.

Dimension lifting is tested in:

- `cfd-1d::geometry_integration_demo` — full pipeline generate→solve→overlay→JSON.
- `cfd-2d::csg_cfd_simulation` — CSG-built channel into 2-D solver, pressure drop
  compared against 1-D network under same geometry.
- `dimension_scenarios_plots` — 1-D/2-D/3-D per-edge pressure and WSS comparison,
  quantifying the $L_2$ deviation as a function of $Re$.

Pressure-drop lifting errors originate from (a) entrance-length effects missed by
1-D fully-developed assumption, (b) separation at sharp turns (2-D and 3-D
produce recirculation not modeled in 1-D), and (c) stenosis-induced jet with
vena contracta contraction — diagnosed by plotting velocity profiles at
discrete cross-sections.

### API — Full Pipeline

```rust
pub trait GeometryPipeline<F: FloatElement> {
    fn build_mesh(&self, primitive: &Mesh, density: &DensitySpec) -> Result<Mesh, CfdError>;
    fn wire_boundaries(&self, mesh: &Mesh) -> Result<Vec<Boundary>, CfdError>;
    fn submit(&self, mesh: &Mesh, boundaries: &[Boundary]) -> Result<SolverHandle, CfdError>;
}
```

CFDrs hides each step behind this trait so that the swap is local: replacing
the mesh backend (from `csgrs` to `hephaestus` GPU) is a single `build_mesh`
implementation change.

---

## 8 — CSG-CFD Coupling — End-to-End

The full pipeline from geometry to posterior is:

```
  CSG primitive  -->  CSG bool op  -->  Mesh  -->  Boundary wiring  -->  Solver  -->  Posterior
     (SDF)             (min/max/-)     (MC / Tet)   (Interface)          (NS)       (WSS, σ, H)
```

CFDrs's design guarantees that geometry work precedes any solver dispatch.
Stenosis screening (`cfd-1d::screening::shear_violations`) walks a chip graph
whose nodes are `NetworkBlueprint` entries (each a venturi / bifurcation /
straight segment built from CSG primitives) and aggregates hemolysis and shear
violations branch-wise.

---

## Examples Referenced by This Chapter

- [Example: csg_primitives_demo](examples/csg_primitives_demo.md) — primitive catalogue walkthrough (box, cylinder, sphere-sections, torus).
- [Example: csg_operations](examples/csg_operations.md) — union / intersection / difference over CSG primitives; SDF evaluation at probe points.
- [Example: csg_cfd_simulation](examples/csg_cfd_simulation.md) — CSG-built channel integrated into the 2-D structured solver with pressure-drop validation.
- [Example: mesh_3d_integration](examples/mesh_3d_integration.md) — 3-D tetrahedral mesh from a CSG boundary; quality metric reporting.
- [Example: dimension_scenarios_plots](examples/dimension_scenarios_plots.md) — 1-D / 2-D / 3-D fidelity ladder for a single venturi geometry.

## Further Reading

- `cfd-schematics` source: `crates/cfd-schematics/src/{primitives,sdf,operations,quality}.rs`.
- [Biomedical flows](biomedical_flows.md) for the geometries that ship pre-built.
- [Turbulence, Multiphase, and Cavitation](turbulence_multiphase.md) for the boundary conditions that the mesh serves.
- [Foundations](foundations.md) for the Atlas geometry type map.
- Hart, "Sphere tracing: a geometric method for the antialiased ray tracing of implicit surfaces", Visual Computer 12, 527 (1996) — sphere tracing for SDF queries.
- Lorensen & Cline, SIGGRAPH 21(4), 163 (1987) — marching cubes.
- Ju et al., SIGGRAPH 21(3), 201 (2002) — dual contouring.
- P. Frey & George, *Mesh Generation*, 2nd ed., Wiley (2008) — Delaunay refinement and quality theory.
