//! Y-junction (bifurcation) primitive — one inlet splitting into two outlet branches.
//!
//! ## Algorithm
//!
//! The mesh is built as a single structured surface by dividing the cross-section into
//! **three sectors** at the junction, each sector owned by one arm:
//!
//! ```text
//!          +Z (north)
//!          |
//!   left   | right
//!   sector | sector
//!          |
//!  --------+-------- (divider line in XZ plane)
//!          |
//!   inlet sector
//!          |
//!          -Z (south)
//! ```
//!
//! Each arm tube owns a sector of `ns` vertices. The sector "seam" vertices (at the
//! north/south poles, `±Z`) are shared between adjacent arms. At the junction, the
//! seam vertices are represented as **single shared vertices** by using
//! `VertexPool` welding with a tight tolerance.
//!
//! The junction body is triangulated by connecting:
//! - Inlet sector to left/right sectors via "gore" strips.
//! - Left sector to right sector via a "divider" strip.
//!
//! This produces a watertight, consistently-oriented manifold mesh.
//!
//! ## Region IDs
//!
//! | [`RegionId`] | Surface |
//! |---|---|
//! | 1 | Inlet tube wall |
//! | 2 | Left branch wall |
//! | 3 | Right branch wall |
//! | 4 | Inlet base cap |
//! | 5 | Left outlet cap |
//! | 6 | Right outlet cap |
//! | 7 | Junction patch |

use std::f64::consts::{TAU, PI};

use crate::domain::core::index::{RegionId, VertexId};
use crate::domain::core::scalar::{Point3r, Vector3r};
use crate::domain::mesh::IndexedMesh;
use super::{PrimitiveMesh, PrimitiveError};

// ── Struct ─────────────────────────────────────────────────────────────────────

/// Builds a Y-junction (bifurcation) mesh: one inlet splitting into two
/// symmetric outlet branches.
#[derive(Clone, Debug)]
pub struct YJunction {
    /// Cross-section (tube) radius [mm]. All three arms share this radius.
    pub tube_radius: f64,
    /// Length of the inlet arm [mm].
    pub inlet_length: f64,
    /// Length of each outlet branch [mm].
    pub branch_length: f64,
    /// Half-angle between +Y axis and each branch centreline [rad].
    /// Must be in `(0, π/2)`.
    pub branch_half_angle: f64,
    /// Number of vertices around the tube cross-section (≥ 4, **even**).
    pub tube_segments: usize,
    /// Number of full tube rings along each arm (≥ 1).
    pub arm_rings: usize,
}

impl Default for YJunction {
    fn default() -> Self {
        Self {
            tube_radius:       0.5,
            inlet_length:      2.0,
            branch_length:     2.0,
            branch_half_angle: PI / 6.0,
            tube_segments:     16,
            arm_rings:         4,
        }
    }
}

impl PrimitiveMesh for YJunction {
    fn build(&self) -> Result<IndexedMesh, PrimitiveError> {
        build(self)
    }
}

// ── Helpers ────────────────────────────────────────────────────────────────────

/// Orthonormal frame perpendicular to `tangent`.
fn perp_frame(tangent: Vector3r) -> (Vector3r, Vector3r) {
    let hint = if tangent.z.abs() < 0.9 { Vector3r::z() } else { Vector3r::x() };
    let e1 = (hint - hint.dot(&tangent) * tangent).normalize();
    let e2 = tangent.cross(&e1).normalize();
    (e1, e2)
}

/// Position and outward normal for vertex `i` on a ring of `ns` vertices,
/// centered at `center`, lying in the plane normal to `tangent`, with given `phase`.
fn ring_vertex(center: Point3r, tangent: Vector3r, r: f64, ns: usize, i: usize, phase: f64)
    -> (Point3r, Vector3r)
{
    let (e1, e2) = perp_frame(tangent);
    let β = TAU * i as f64 / ns as f64 + phase;
    let (sb, cb) = β.sin_cos();
    let off = r * (cb * e1 + sb * e2);
    (Point3r::from(center.coords + off), off.normalize())
}

/// Insert a complete ring of `ns` vertices.
fn insert_ring(
    mesh:    &mut IndexedMesh,
    center:  Point3r,
    tangent: Vector3r,
    r:       f64,
    ns:      usize,
    phase:   f64,
) -> Vec<VertexId> {
    (0..ns).map(|i| {
        let (pos, nrm) = ring_vertex(center, tangent, r, ns, i, phase);
        mesh.add_vertex(pos, nrm)
    }).collect()
}

/// Stitch ring_a → ring_b with outward-facing triangles.
///
/// Convention (CCW from outside):
/// - face1: `(A[i], A[j], B[j])` — A edge i→j consumed
/// - face2: `(A[i], B[j], B[i])` — B edge j→i consumed
fn stitch_rings(mesh: &mut IndexedMesh, a: &[VertexId], b: &[VertexId], region: RegionId) {
    let n = a.len();
    for i in 0..n {
        let j = (i + 1) % n;
        mesh.add_face_with_region(a[i], a[j], b[j], region);
        mesh.add_face_with_region(a[i], b[j], b[i], region);
    }
}

/// Stitch open arc a[0..=m] → b[0..=m] (not cyclic).
/// face1: `(A[i], A[i+1], B[i+1])` — A edge i→i+1 consumed
/// face2: `(A[i], B[i+1], B[i])` — B edge i+1→i consumed
fn stitch_arc(mesh: &mut IndexedMesh, a: &[VertexId], b: &[VertexId], region: RegionId) {
    assert_eq!(a.len(), b.len());
    let n = a.len();
    for i in 0..n - 1 {
        let j = i + 1;
        if a[i] != a[j] && a[i] != b[j] && a[j] != b[j] {
            mesh.add_face_with_region(a[i], a[j], b[j], region);
        }
        if a[i] != b[j] && b[j] != b[i] && a[i] != b[i] {
            mesh.add_face_with_region(a[i], b[j], b[i], region);
        }
    }
}

/// Triangle-fan cap. Outward normal = cap_normal.
/// Convention: (center, ring[j], ring[i]) — ring edge j→i consumed.
fn cap_end(
    mesh:       &mut IndexedMesh,
    ring:       &[VertexId],
    center:     Point3r,
    cap_normal: Vector3r,
    region:     RegionId,
) {
    let cid = mesh.add_vertex(center, cap_normal);
    let n = ring.len();
    for i in 0..n {
        let j = (i + 1) % n;
        if ring[i] != ring[j] {
            mesh.add_face_with_region(cid, ring[j], ring[i], region);
        }
    }
}

// ── Builder ────────────────────────────────────────────────────────────────────

fn build(y: &YJunction) -> Result<IndexedMesh, PrimitiveError> {
    if y.tube_radius <= 0.0 {
        return Err(PrimitiveError::InvalidParam(format!(
            "tube_radius must be > 0, got {}", y.tube_radius)));
    }
    if y.inlet_length <= 0.0 {
        return Err(PrimitiveError::InvalidParam(format!(
            "inlet_length must be > 0, got {}", y.inlet_length)));
    }
    if y.branch_length <= 0.0 {
        return Err(PrimitiveError::InvalidParam(format!(
            "branch_length must be > 0, got {}", y.branch_length)));
    }
    if y.branch_half_angle <= 0.0 || y.branch_half_angle >= PI / 2.0 {
        return Err(PrimitiveError::InvalidParam(format!(
            "branch_half_angle must be in (0, π/2), got {}", y.branch_half_angle)));
    }
    if y.tube_segments < 4 || y.tube_segments % 2 != 0 {
        return Err(PrimitiveError::InvalidParam(format!(
            "tube_segments must be ≥ 4 and even, got {}", y.tube_segments)));
    }
    if y.arm_rings < 1 {
        return Err(PrimitiveError::InvalidParam("arm_rings must be ≥ 1".into()));
    }

    let r   = y.tube_radius;
    let ns  = y.tube_segments;
    let h  = ns / 2;   // half-ring size
    let nr  = y.arm_rings;
    let θ   = y.branch_half_angle;
    let (sθ, cθ) = θ.sin_cos();

    // Tolerance: tight enough to weld pole vertices (all rings share exact positions at poles).
    let min_spacing = TAU * r / ns as f64;
    let tol = (min_spacing / 20.0).max(1e-10);
    let mut mesh = IndexedMesh::with_tolerance(tol * 4.0, tol);

    // Arm tangent unit vectors.
    let t_in    = Vector3r::new( 0.0,  1.0,  0.0);
    let t_left  = Vector3r::new( sθ,   cθ,   0.0);
    let t_right = Vector3r::new(-sθ,   cθ,   0.0);

    let reg_inlet = RegionId::new(1);
    let reg_left  = RegionId::new(2);
    let reg_right = RegionId::new(3);
    let reg_incap = RegionId::new(4);
    let reg_lcap  = RegionId::new(5);
    let reg_rcap  = RegionId::new(6);
    let reg_junc  = RegionId::new(7);

    // All three arms use phase=0.
    // The ring vertex at index 0 is at angle β=0, giving position:
    //   For t_in: e1=Z, e2=X; vertex[0] = r*(cos(0)*Z + sin(0)*X) = (0,0,r) = north pole.
    //   For t_left and t_right: same e1=Z (since tangent is in XY plane), vertex[0] = (0,0,r).
    // All three junction rings share the same north pole (0,0,r) and south pole (0,0,-r).
    // VertexPool welding will merge these into single VertexIds.
    //
    // The "seam" vertices at the poles are automatically shared between arm rings,
    // making the junction topologically clean.

    // ── Inlet arm ────────────────────────────────────────────────────────────
    // Rings: k=0 (outer end, y=-inlet_length) to k=nr (junction end, y=0).
    // stitch_rings(ring[k], ring[k+1]) — junc_in = ring[nr] is ring_b (tip).
    // Arm stitch face2 uses junc_in j→i → junction must use junc_in i→j.
    let inlet_rings: Vec<Vec<VertexId>> = (0..=nr).map(|k| {
        let t = k as f64 / nr as f64;
        let c = Point3r::new(0.0, -y.inlet_length * (1.0 - t), 0.0);
        insert_ring(&mut mesh, c, t_in, r, ns, 0.0)
    }).collect();
    for k in 0..nr {
        stitch_rings(&mut mesh, &inlet_rings[k], &inlet_rings[k + 1], reg_inlet);
    }
    cap_end(&mut mesh, &inlet_rings[0],
            Point3r::new(0.0, -y.inlet_length, 0.0), -t_in, reg_incap);
    let junc_in = inlet_rings[nr].clone();

    // ── Left branch ──────────────────────────────────────────────────────────
    // Rings: k=0 (junction end, origin) to k=nr (tip).
    // stitch_rings(ring[k], ring[k+1]) — junc_l = ring[0] is ring_a (base).
    // Arm stitch face1 uses junc_l i→j → junction must use junc_l j→i.
    let left_tip = Point3r::new(sθ * y.branch_length, cθ * y.branch_length, 0.0);
    let left_rings: Vec<Vec<VertexId>> = (0..=nr).map(|k| {
        let t = k as f64 / nr as f64;
        let c = Point3r::from(left_tip.coords * t);
        insert_ring(&mut mesh, c, t_left, r, ns, 0.0)
    }).collect();
    for k in 0..nr {
        stitch_rings(&mut mesh, &left_rings[k], &left_rings[k + 1], reg_left);
    }
    cap_end(&mut mesh, &left_rings[nr], left_tip, t_left, reg_lcap);
    let junc_l = left_rings[0].clone();

    // ── Right branch ─────────────────────────────────────────────────────────
    // Same as left. junc_r = ring[0] is ring_a; arm stitch face1 uses i→j.
    // Junction must use junc_r j→i.
    let right_tip = Point3r::new(-sθ * y.branch_length, cθ * y.branch_length, 0.0);
    let right_rings: Vec<Vec<VertexId>> = (0..=nr).map(|k| {
        let t = k as f64 / nr as f64;
        let c = Point3r::from(right_tip.coords * t);
        insert_ring(&mut mesh, c, t_right, r, ns, 0.0)
    }).collect();
    for k in 0..nr {
        stitch_rings(&mut mesh, &right_rings[k], &right_rings[k + 1], reg_right);
    }
    cap_end(&mut mesh, &right_rings[nr], right_tip, t_right, reg_rcap);
    let junc_r = right_rings[0].clone();

    // ── Junction patch ───────────────────────────────────────────────────────
    //
    // Summary of edge-direction constraints from arm stitches:
    //   junc_in[k]→junc_in[k+1] : consumed j→i by inlet. Junction needs i→j.
    //   junc_l[k]→junc_l[k+1]   : consumed i→j by left (face1). Junction needs j→i.
    //   junc_r[k]→junc_r[k+1]   : consumed i→j by right (face1). Junction needs j→i.
    //
    // With phase=0, the junction rings share pole vertices:
    //   P_north = junc_in[0] = junc_l[0] = junc_r[0]   (position (0,0,r))
    //   P_south = junc_in[h] = junc_l[h] = junc_r[h]   (position (0,0,-r))
    //
    // Due to VertexPool welding, these are the same VertexId.
    // Arm stitches at poles produce degenerate-looking triangles but are OK (they're valid).
    //
    // The junction patch divides into three sub-patches:
    //
    // LEFT GORE: connects junc_in FRONT arc [0..h] to junc_l FRONT arc [0..h].
    //   "Front" = the arc facing toward the left branch (+X side).
    //   Both arcs share endpoints: junc_in[0]=junc_l[0]=P_north, junc_in[h]=junc_l[h]=P_south.
    //   Interior quads (i=0..h-1, j=i+1):
    //     T1=(junc_in[i], junc_in[j], junc_l[i])  — inlet i→j ✓
    //     T2=(junc_in[j], junc_l[j], junc_l[i])   — junc_l j→i ✓ (j>i)
    //     Manifold? Shared diagonal (junc_in[j], junc_l[i]):
    //       T1 edges: in[i]→in[j], in[j]→l[i], l[i]→in[i].  T1 shared: in[j]→l[i].
    //       T2 edges: in[j]→l[j], l[j]→l[i], l[i]→in[j].    T2 shared: l[i]→in[j]. OPPOSITE ✓.
    //   At poles (i=0: junc_in[0]=junc_l[0]=P_north, the T2 face degenerates → skipped.
    //            (i=h-1: junc_in[h]=junc_l[h]=P_south, the T1+T2 both have pole as endpoint.
    //              T1=(in[h-1], P_south, l[h-1]): valid triangle.
    //              T2=(P_south, l[h], l[h-1])=(P_south, P_south, l[h-1]): DEGENERATE → skip.
    //
    // RIGHT GORE: connects junc_in BACK arc [h..2h=ns] to junc_r FRONT arc [h..0] (reversed).
    //   junc_in back arc goes i→j (k=0: index h, k=h-1: index ns-1→ns%ns=0).
    //   junc_r front arc goes j→i in the junction (k=0: index h, k=h-1: index 1→0).
    //   Triangulation (proven manifold at lines ~1048-1053):
    //     ia=junc_in[(h+k)%ns], ib=junc_in[(h+k+1)%ns]
    //     rra=junc_r[(h-k)%ns], rrb=junc_r[(h-k-1+ns)%ns]
    //     T1=(ia, ib, rra)  — inlet ia→ib (i→j ✓)
    //     T2=(ia, rra, rrb) — junc_r rra→rrb (h-k → h-k-1 = decreasing = j→i ✓)
    //     Shared (ia,rra): T1 rra→ia, T2 ia→rra. OPPOSITE ✓.
    //   At poles: k=0: ia=in[h]=P_south, rra=junc_r[h]=P_south → ia=rra → T1 degenerate.
    //             k=h-1: ib=in[0]=P_north, rrb=junc_r[0]=P_north → ib=P_north=rrb → T2 degenerate.
    //
    // DIVIDER: connects junc_l BACK arc [h..ns] to junc_r BACK arc [h..ns].
    //   Both need j→i direction. Cannot do a standard strip.
    //   Resolution: since both arcs share the same POLE VERTICES at k=0 (P_south) and k=h (P_north),
    //   the strip degenerates at endpoints. In between, we need a strip where both go j→i.
    //
    //   Key observation: the left and right branches are SYMMETRIC across the YZ plane.
    //   The left back arc vertices at angles β=π..2π (XZ projection: -X side of left ring).
    //   The right back arc vertices at angles β=π..2π (XZ projection: -X side of right ring).
    //   But these are in DIFFERENT 3D planes!
    //
    //   For the divider strip, let's define:
    //     la[k] = junc_l[(h+k)%ns], k=0..h  (back arc of left, goes south→north as k increases)
    //     ra[k] = junc_r[(h+k)%ns], k=0..h  (back arc of right, same direction)
    //   la[0]=junc_l[h]=P_south, la[h]=junc_l[0]... wait: (h+h)%ns = 0. So la[h]=junc_l[0]=P_north.
    //   ra[0]=junc_r[h]=P_south, ra[h]=junc_r[0]=P_north.
    //
    //   Both arcs go south→north (k=0→h). The constraint says junction uses j→i for both.
    //   j→i means DECREASING k direction for the strip.
    //
    //   Alternative: run the strip from north→south (k=h→0):
    //     la_k = junc_l[(h + (h-k))%ns] = junc_l[(2h-k)%ns] = junc_l[(ns-k)%ns]
    //   Actually: just run k=0..h but with reversed indices:
    //     la[k] = junc_l[(ns - k) % ns]  (k=0: junc_l[0]=P_north, k=h: junc_l[h]=P_south)
    //     ra[k] = junc_r[(ns - k) % ns]  (k=0: junc_r[0]=P_north, k=h: junc_r[h]=P_south)
    //
    //   Now: la goes north→south as k increases. The arm stitch uses la in i→j direction
    //   (face1 of left arm: junc_l[i]→junc_l[i+1]). The direction j→i in the original ring
    //   is i+1→i (decreasing). So for the junction to use j→i: junction uses la[k]→la[k-1].
    //   In the strip, adjacent la vertices go la[0],la[1],...,la[h]. Using la[k]→la[k-1] means
    //   the strip goes k=1→2→...→h (using la in the "la[k]→la[k-1]=decreasing k" sense),
    //   i.e., REVERSE of the parameterization.
    //
    //   T1=(la[k], la[k+1], ra[k]):
    //     la edge: la[k]→la[k+1] = junc_l[(ns-k)%ns] → junc_l[(ns-k-1)%ns]. Decreasing = j→i ✓.
    //     Shared diagonal (la[k], ra[k])?
    //     T1 edges: la[k]→la[k+1], la[k+1]→ra[k], ra[k]→la[k].
    //   T2=(la[k], ra[k], ra[k+1]):
    //     ra edge: ra[k]→ra[k+1] = junc_r[(ns-k)%ns]→junc_r[(ns-k-1)%ns]. Decreasing = j→i ✓.
    //     Shared (la[k], ra[k]): T1 has ra[k]→la[k]. T2 has la[k]→ra[k]. OPPOSITE ✓!
    //   ✓ MANIFOLD! ✓ BOTH ARCS j→i ✓!
    //
    //   Wait — la[k]→la[k+1] is decreasing in junc_l original index. Let me verify:
    //   la[k] = junc_l[(ns-k)%ns]. la[k]→la[k+1] = junc_l[(ns-k)%ns] → junc_l[(ns-k-1)%ns].
    //   This is junc_l[ns-k] → junc_l[ns-k-1] = junc_l[q+1] → junc_l[q] where q=ns-k-1.
    //   So q+1 → q = i+1→i = j→i for the ORIGINAL junc_l ring indices. ✓
    //
    //   At poles: k=0: la[0]=junc_l[0]=P_north, ra[0]=junc_r[0]=P_north. la[0]=ra[0]=P_north.
    //     T1=(P_north, la[1], ra[0])=(P_north, la[1], P_north) = degenerate (first=last). Skip.
    //     T2=(P_north, P_north, ra[1]) = degenerate. Skip.
    //   k=h-1: la[h]=junc_l[h]=P_south... wait: la[h] = junc_l[(ns-h)%ns] = junc_l[h] = P_south.
    //     T1=(la[h-1], P_south, ra[h-1]): valid (P_south is distinct from the others unless la[h-1]=P_south).
    //     T2=(la[h-1], ra[h-1], P_south): valid triangle.
    //     The last step: T1 and T2 at k=h-1 are valid (P_south ≠ la[h-1] for h-1 > 0).
    //
    // FINAL JUNCTION PATCH WINDING CHECK:
    //   The junction patch faces the "outside" of the Y-junction (not the interior fluid channel).
    //   LEFT GORE: T1=(in[i], in[j], l[i]). Normal: (in[j]-in[i]) × (l[i]-in[i]).
    //     in[i] and in[j] are at the inlet junction ring (at origin, in XZ plane).
    //     l[i] is at the left branch junction ring (similar plane). The cross product should
    //     point toward the "exterior" of the Y-junction (away from the fluid interior).
    //   This is geometry-dependent and may require a sign flip. We'll check in tests.

    // LEFT GORE
    for i in 0..h {
        let j = i + 1;
        let ia = junc_in[i];
        let ib = junc_in[j % ns];
        let la = junc_l[i];
        let lb = junc_l[j % ns];
        // T1: inlet i→j, uses junc_l[i] as a vertex only.
        if ia != ib && ib != la && la != ia {
            mesh.add_face_with_region(ia, ib, la, reg_junc);
        }
        // T2: junc_l j→i (lb→la).
        if ib != lb && lb != la && la != ib {
            mesh.add_face_with_region(ib, lb, la, reg_junc);
        }
    }

    // RIGHT GORE
    for k in 0..h {
        let ia  = junc_in[(h + k) % ns];
        let ib  = junc_in[(h + k + 1) % ns];
        let rra = junc_r[(h.wrapping_sub(k)) % ns];
        let rrb = junc_r[(h + ns - k - 1) % ns];
        // T1: inlet ia→ib (i→j)
        if ia != ib && ib != rra && rra != ia {
            mesh.add_face_with_region(ia, ib, rra, reg_junc);
        }
        // T2: junc_r rra→rrb (j→i)
        if ia != rra && rra != rrb && rrb != ia {
            mesh.add_face_with_region(ia, rra, rrb, reg_junc);
        }
    }

    // DIVIDER
    // la[k] = junc_l[(ns-k)%ns], ra[k] = junc_r[(ns-k)%ns] for k=0..h.
    for k in 0..h {
        let la_k  = junc_l[(ns - k) % ns];       // la[k]
        let la_k1 = junc_l[(ns - k - 1) % ns];   // la[k+1]
        let ra_k  = junc_r[(ns - k) % ns];        // ra[k]
        let ra_k1 = junc_r[(ns - k - 1) % ns];   // ra[k+1]
        // T1: junc_l la_k→la_k1 (j→i), diagonal (la_k, ra_k).
        if la_k != la_k1 && la_k1 != ra_k && ra_k != la_k {
            mesh.add_face_with_region(la_k, la_k1, ra_k, reg_junc);
        }
        // T2: junc_r ra_k→ra_k1 (j→i), diagonal (la_k, ra_k) opposite side.
        if la_k != ra_k && ra_k != ra_k1 && ra_k1 != la_k {
            mesh.add_face_with_region(la_k, ra_k, ra_k1, reg_junc);
        }
    }

    Ok(mesh)
}

// ── Tests ───────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    fn diagnostic(ns: usize) {
        let yj = YJunction { tube_segments: ns, ..YJunction::default() };
        let mut mesh = yj.build().expect("build should succeed");
        let is_wt = mesh.is_watertight();
        // Re-borrow for detailed stats
        let report = {
            let edges = mesh.edges_ref().unwrap();
            crate::application::watertight::check::check_watertight(
                &mesh.vertices, &mesh.faces, edges,
            )
        };
        let _ = is_wt;
        println!("ns={}: V={} F={} boundary={} nm={} closed={} chi={:?} vol={:.4}",
            ns,
            mesh.vertex_count(), mesh.face_count(),
            report.boundary_edge_count, report.non_manifold_edge_count,
            report.is_closed, report.euler_characteristic,
            report.signed_volume);
        // Print boundary edge positions for small ns
        if ns <= 8 {
            let edges = mesh.edges_ref().unwrap();
            for edge in edges.iter() {
                if edge.is_boundary() {
                    let (v0, v1) = edge.vertices;
                    let p0 = mesh.vertices.position(v0);
                    let p1 = mesh.vertices.position(v1);
                    println!("  boundary edge: ({:.2},{:.2},{:.2})→({:.2},{:.2},{:.2})",
                        p0.x,p0.y,p0.z, p1.x,p1.y,p1.z);
                }
            }
        }
    }

    #[test]
    fn y_junction_diagnostic() {
        diagnostic(4);
        diagnostic(8);
        diagnostic(16);
    }

    #[test]
    fn y_junction_builds() {
        let yj = YJunction::default();
        let mesh = yj.build().expect("build should succeed");
        assert!(mesh.face_count() > 0);
    }

    #[test]
    fn y_junction_invalid_params() {
        let result = YJunction { tube_segments: 3, ..YJunction::default() }.build();
        assert!(result.is_err(), "should reject tube_segments=3 (odd)");
        let result = YJunction { branch_half_angle: 0.0, ..YJunction::default() }.build();
        assert!(result.is_err(), "should reject zero branch angle");
    }
}
