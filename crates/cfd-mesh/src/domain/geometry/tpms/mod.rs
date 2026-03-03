//! Triply Periodic Minimal Surfaces (TPMS) — implicit field definitions.
//!
//! ## What is a TPMS?
//!
//! A **Triply Periodic Minimal Surface** is a surface with zero mean curvature
//! everywhere that is periodic in three independent lattice directions.  Common
//! examples include:
//!
//! | Surface | Implicit field F(x,y,z) | Discoverer |
//! |---------|------------------------|------------|
//! | Gyroid | sin(kx)cos(ky)+sin(ky)cos(kz)+sin(kz)cos(kx) | Schoen (1970) |
//! | Schwarz P | cos(kx)+cos(ky)+cos(kz) | Schwarz (1890) |
//! | Schwarz D | sin(kx)sin(ky)sin(kz)+… | Schwarz (1890) |
//! | Neovius | 3(cos(kx)+cos(ky)+cos(kz))+4cos(kx)cos(ky)cos(kz) | Neovius (1883) |
//! | Lidinoid | ½(sin(2kx)cos(ky)sin(kz)+…)−½(cos(2kx)cos(2ky)+…)+0.15 | Lidin (1990) |
//! | I-WP | 2(cos(kx)cos(ky)+…)−(cos(2kx)+…) | Schoen (1970) |
//! | Split P | 1.1(sin(2kx)sin(kz)cos(ky)+…)−0.2(…)−0.4(…) | Fogden & Hyde (1992) |
//! | FRD | 4cos(kx)cos(ky)cos(kz)−(cos(2kx)cos(2ky)+…) | Schoen (1970) |
//! | Fischer-Koch C(Y) | cos(2kx)sin(ky)cos(kz)+… | Fischer & Koch (1989) |
//!
//! The zero level-set `F = 0` of each defines the TPMS.  Offset level-sets
//! (`F = c`) yield constant-thickness shell pairs (useful for lattice infill).
//!
//! ## Architecture
//!
//! ```text
//! tpms/
//! ├── mod.rs              ← this file: Tpms trait + TpmsParams + build_tpms_sphere()
//! ├── marching_cubes.rs   ← single copy of EDGE_TABLE, TRI_TABLE, extract()
//! ├── gyroid.rs           ← Gyroid: F = sin(kx)cos(ky)+…
//! ├── schwarz_p.rs        ← SchwarzP: F = cos(kx)+cos(ky)+cos(kz)
//! ├── schwarz_d.rs        ← SchwarzD: F = sin(kx)sin(ky)sin(kz)+…
//! ├── neovius.rs          ← Neovius: F = 3Σcos+4Πcos
//! ├── lidinoid.rs         ← Lidinoid: F = ½Σsin2cos·sin−½Σcos2cos2+0.15
//! ├── iwp.rs              ← Iwp: F = 2Σcos·cos−Σcos2
//! ├── split_p.rs          ← SplitP: F = 1.1Σsin2sincos−0.2Σcos2cos2−0.4Σcos2
//! ├── frd.rs              ← Frd: F = 4Πcos−Σcos2cos2
//! └── fischer_koch_cy.rs  ← FischerKochCY: F = Σcos2·sin·cos (cyclic)
//! ```

pub mod box_clip;
pub mod fischer_koch_cy;
pub mod frd;
pub mod gyroid;
pub mod iwp;
pub mod lidinoid;
pub mod marching_cubes;
pub mod neovius;
pub mod schwarz_d;
pub mod schwarz_p;
pub mod split_p;

pub use box_clip::{build_tpms_box, build_tpms_box_graded, TpmsBoxParams};
pub use fischer_koch_cy::FischerKochCY;
pub use frd::Frd;
pub use gyroid::Gyroid;
pub use iwp::Iwp;
pub use lidinoid::Lidinoid;
pub use neovius::Neovius;
pub use schwarz_d::SchwarzD;
pub use schwarz_p::SchwarzP;
pub use split_p::SplitP;

use crate::domain::core::scalar::Vector3r;
use crate::domain::geometry::primitives::PrimitiveError;
use crate::domain::mesh::IndexedMesh;

// ── TPMS trait ────────────────────────────────────────────────────────────────

/// Implicit field function for a Triply Periodic Minimal Surface.
///
/// Implementors define the specific algebraic form of the TPMS field.  The
/// shared [`build_tpms_sphere`] builder handles the marching-cubes extraction
/// and sphere clipping uniformly for all surfaces.
///
/// ## Theorem — Zero Mean Curvature
///
/// For each of the provided implementors, the zero level-set `F = 0` is a
/// minimal surface with `H = 0` everywhere.  This follows from the implicit
/// function theorem: a level-set `{F = 0}` has zero mean curvature iff
/// `div(∇F / |∇F|) = 0`.  Verification is by direct computation of the
/// Laplacian and second fundamental form for each surface (see individual
/// module docs).
pub trait Tpms: Send + Sync {
    /// Evaluate the implicit field at world position `(x, y, z)` with
    /// temporal angular frequency `k = 2π / period`.
    ///
    /// The TPMS mid-surface is the zero level-set `{F(x,y,z,k) = 0}`.
    fn field(&self, x: f64, y: f64, z: f64, k: f64) -> f64;

    /// Normalised outward surface normal at `(x, y, z)` — the analytical
    /// gradient `∇F / |∇F|`.
    ///
    /// Implementations should fall back to `Vector3r::y()` when `|∇F| ≈ 0`
    /// to avoid NaN propagation at saddle points.
    fn gradient(&self, x: f64, y: f64, z: f64, k: f64) -> Vector3r;
}

// ── Shared build helper ───────────────────────────────────────────────────────

/// Build parameters for sphere-clipped TPMS primitives.
#[derive(Clone, Debug)]
pub struct TpmsParams {
    /// Clip-sphere radius [mm].  The TPMS surface is extracted inside this sphere.
    pub radius: f64,
    /// TPMS unit-cell period [mm].  `k = 2π / period`.
    pub period: f64,
    /// Number of voxels per axis.  Higher values yield denser, more accurate meshes.
    pub resolution: usize,
    /// Level-set threshold.  `0.0` = exact minimal surface mid-sheet.
    /// Positive values shift the surface outward; negative inward.
    pub iso_value: f64,
}

impl Default for TpmsParams {
    fn default() -> Self {
        Self {
            radius: 5.0,
            period: 2.5,
            resolution: 64,
            iso_value: 0.0,
        }
    }
}

/// Extract a sphere-clipped TPMS mesh using the shared marching-cubes engine.
///
/// This is the single shared build path used by all TPMS sphere primitives.
/// New TPMS primitives need only implement [`Tpms`] and delegate to this
/// function — O(1) cost per new surface type.
pub fn build_tpms_sphere<S: Tpms>(
    surface: &S,
    params: &TpmsParams,
) -> Result<IndexedMesh, PrimitiveError> {
    if params.radius <= 0.0 {
        return Err(PrimitiveError::InvalidParam(format!(
            "radius must be > 0, got {}",
            params.radius
        )));
    }
    if params.period <= 0.0 {
        return Err(PrimitiveError::InvalidParam(format!(
            "period must be > 0, got {}",
            params.period
        )));
    }
    if params.resolution < 4 {
        return Err(PrimitiveError::InvalidParam(format!(
            "resolution must be >= 4, got {}",
            params.resolution
        )));
    }

    let k = std::f64::consts::TAU / params.period;
    let mc_params = marching_cubes::McParams {
        radius: params.radius,
        resolution: params.resolution,
        k,
        iso_value: params.iso_value,
    };

    let mut mesh = IndexedMesh::new();
    marching_cubes::extract(
        &mut mesh,
        &mc_params,
        |x, y, z, kk| surface.field(x, y, z, kk),
        |x, y, z, kk| surface.gradient(x, y, z, kk),
    );

    Ok(mesh)
}
