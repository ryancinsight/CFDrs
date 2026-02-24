//! Triangle quality measurement.
//!
//! # Radius-Edge Ratio
//!
//! The primary quality metric for Ruppert's algorithm is the radius-edge ratio:
//!
//! $$\rho(t) = \frac{R(t)}{l_{\min}(t)}$$
//!
//! where $R(t)$ is the circumradius and $l_{\min}(t)$ is the shortest edge
//! length of triangle $t$.
//!
//! # Theorem — Minimum Angle Relationship
//!
//! **Statement**: For a triangle $t$,
//!
//! $$\rho(t) = \frac{1}{2 \sin(\alpha_{\min})}$$
//!
//! where $\alpha_{\min}$ is the minimum interior angle.
//!
//! **Proof sketch**: By the law of sines, each edge $a$ opposite angle $A$
//! satisfies $\frac{a}{\sin A} = 2R$.  The shortest edge is opposite the
//! smallest angle, so $l_{\min} = 2R \sin(\alpha_{\min})$.
//! Hence $\rho = R / (2R \sin \alpha_{\min}) = 1 / (2 \sin \alpha_{\min})$.
//!
//! # Consequence
//!
//! A bound $\rho \leq B$ implies $\alpha_{\min} \geq \arcsin(1/(2B))$.
//!
//! | $B$ | $\alpha_{\min}$ (degrees) |
//! |-----|---------------------------|
//! | $\sqrt{2} \approx 1.414$ | 20.7° |
//! | $1.0$ | 30.0° |
//! | $0.58$ | 59.3° (near equilateral) |

use crate::application::delaunay::pslg::vertex::PslgVertex;
use crate::domain::core::scalar::Real;

/// Quality metrics for a single triangle.
#[derive(Debug, Clone, Copy)]
pub struct TriangleQuality {
    /// The circumradius.
    pub circumradius: Real,
    /// The shortest edge length.
    pub shortest_edge: Real,
    /// The radius-edge ratio $R / l_{\min}$.
    pub radius_edge_ratio: Real,
    /// The minimum interior angle in radians.
    pub min_angle_rad: Real,
    /// The area.
    pub area: Real,
}

impl TriangleQuality {
    /// Compute quality metrics for the triangle `(a, b, c)`.
    pub fn compute(a: &PslgVertex, b: &PslgVertex, c: &PslgVertex) -> Self {
        let (abx, aby) = (b.x - a.x, b.y - a.y);
        let (bcx, bcy) = (c.x - b.x, c.y - b.y);
        let (cax, cay) = (a.x - c.x, a.y - c.y);

        let lab = (abx * abx + aby * aby).sqrt();
        let lbc = (bcx * bcx + bcy * bcy).sqrt();
        let lca = (cax * cax + cay * cay).sqrt();

        let shortest = lab.min(lbc).min(lca);

        // Signed area (CCW positive).
        let area_2 = abx * (c.y - a.y) - aby * (c.x - a.x);
        let area = area_2.abs() * 0.5;

        // Circumradius: R = (a * b * c) / (4 * area)
        let circumradius = if area > 1e-30 {
            (lab * lbc * lca) / (4.0 * area)
        } else {
            Real::INFINITY
        };

        let ratio = if shortest > 1e-30 {
            circumradius / shortest
        } else {
            Real::INFINITY
        };

        // Minimum angle via law of cosines.
        let min_angle = Self::min_interior_angle(lab, lbc, lca);

        Self {
            circumradius,
            shortest_edge: shortest,
            radius_edge_ratio: ratio,
            min_angle_rad: min_angle,
            area,
        }
    }

    /// Check if the triangle quality satisfies the given radius-edge ratio bound.
    pub fn is_good(&self, max_ratio: Real) -> bool {
        self.radius_edge_ratio <= max_ratio
    }

    /// Minimum angle in degrees.
    pub fn min_angle_deg(&self) -> Real {
        self.min_angle_rad.to_degrees()
    }

    /// Compute the minimum interior angle of a triangle given edge lengths.
    fn min_interior_angle(ab: Real, bc: Real, ca: Real) -> Real {
        // Angle at A (between edges AB and CA): cos A = (AB² + CA² - BC²) / (2·AB·CA)
        let cos_a = (ab * ab + ca * ca - bc * bc) / (2.0 * ab * ca);
        let cos_b = (ab * ab + bc * bc - ca * ca) / (2.0 * ab * bc);
        let cos_c = (bc * bc + ca * ca - ab * ab) / (2.0 * bc * ca);

        let angle_a = cos_a.clamp(-1.0, 1.0).acos();
        let angle_b = cos_b.clamp(-1.0, 1.0).acos();
        let angle_c = cos_c.clamp(-1.0, 1.0).acos();

        angle_a.min(angle_b).min(angle_c)
    }
}
