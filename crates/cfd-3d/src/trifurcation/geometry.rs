//! 3D Trifurcation geometry modeling
//!
//! Models a vessel branching into three daughter vessels.
//! Follows the same principles as Bifurcation3D but with three outlets.
//!
//! # Theorem — Generalised Murray’s Law for N-way Branching
//!
//! For a parent vessel of diameter $D_0$ splitting into $N$ daughters $D_1, \ldots, D_N$,
//! the optimal diameter relationship minimising total transport cost is
//!
//! ```text
//! D₀³ = Σᵢ Dᵢ³
//! ```
//!
//! For a symmetric trifurcation ($D_1 = D_2 = D_3 = D_d$) this reduces to
//! $D_0 = 3^{1/3} D_d \approx 1.442 D_d$.
//!
//! # Theorem — Mass Conservation at Junction
//!
//! For incompressible flow, continuity at the junction requires
//!
//! ```text
//! Q₀ = Q₁ + Q₂ + Q₃
//! ```
//!
//! where $Q_i = u_i A_i$ is the volumetric flow rate in each branch.
//! Combined with Hagen–Poiseuille ($Q \propto D^4 \Delta P / L$), this
//! uniquely determines the pressure drop distribution.

use crate::bifurcation::ConicalTransition;
use cfd_core::conversion::SafeFromF64;
use nalgebra::RealField;
use num_traits::{FromPrimitive, ToPrimitive};
use serde::{Deserialize, Serialize};

/// 3D Trifurcation geometry
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TrifurcationGeometry3D<T: cfd_mesh::domain::core::Scalar + RealField + Copy> {
    /// Parent branch diameter [m]
    pub d_parent: T,
    /// Parent branch length [m]
    pub l_parent: T,

    /// Daughter diameters [m]
    pub d_daughters: [T; 3],
    /// Daughter lengths [m]
    pub l_daughters: [T; 3],

    /// Transition zone length [m]
    pub l_transition: T,
    /// Transition model
    pub transition: ConicalTransition<T>,

    /// Angles between daughters and parent axis [radians]
    /// Typically symmetric, e.g., [30°, 0°, -30°]
    pub branching_angles: [T; 3],
}

impl<
        T: cfd_mesh::domain::core::Scalar
            + RealField
            + Copy
            + FromPrimitive
            + ToPrimitive
            + SafeFromF64,
    > TrifurcationGeometry3D<T>
{
    /// Create a symmetric 3D trifurcation
    pub fn symmetric(
        d_parent: T,
        d_daughter: T,
        l_parent: T,
        l_daughter: T,
        l_transition: T,
        spread_angle: T,
    ) -> Self {
        Self {
            d_parent,
            l_parent,
            d_daughters: [d_daughter, d_daughter, d_daughter],
            l_daughters: [l_daughter, l_daughter, l_daughter],
            l_transition,
            transition: ConicalTransition::SmoothCone {
                length: l_transition,
            },
            branching_angles: [spread_angle, T::zero(), -spread_angle],
        }
    }

    /// Verification based on Generalized Murray's Law: D_0^3 = Σ D_i^3
    pub fn murray_law_deviation(&self) -> T {
        let sum_daughters = num_traits::Float::powi(self.d_daughters[0], 3)
            + num_traits::Float::powi(self.d_daughters[1], 3)
            + num_traits::Float::powi(self.d_daughters[2], 3);
        num_traits::Float::powi(self.d_parent, 3) - sum_daughters
    }
}

impl<
        T: cfd_mesh::domain::core::Scalar
            + nalgebra::RealField
            + Copy
            + num_traits::FromPrimitive
            + num_traits::ToPrimitive
            + cfd_core::conversion::SafeFromF64,
    > cfd_mesh::application::delaunay::dim3::sdf::Sdf3D<T> for TrifurcationGeometry3D<T>
{
    fn eval(&self, p: &nalgebra::Point3<T>) -> T {
        use cfd_mesh::application::delaunay::dim3::sdf::FiniteCylinderSdf;
        use nalgebra::{Point3, Vector3};
        
        // Start with parent cylinder distance
        let half = <T as cfd_mesh::domain::core::Scalar>::from_f64(2.0);
        
        let parent = FiniteCylinderSdf::new(
            Point3::new(-self.l_parent, T::zero(), T::zero()),
            Point3::new(T::zero(), T::zero(), T::zero()),
            self.d_parent / half,
        );
        let mut min_val = parent.eval(p);
        
        for i in 0..3 {
            let angle = self.branching_angles[i];
            let dir = Vector3::new(num_traits::Float::cos(angle), num_traits::Float::sin(angle), T::zero());
            let daughter = FiniteCylinderSdf::new(
                Point3::new(T::zero(), T::zero(), T::zero()),
                Point3::new(T::zero(), T::zero(), T::zero()) + dir * self.l_daughters[i],
                self.d_daughters[i] / half,
            );
            min_val = num_traits::Float::min(min_val, daughter.eval(p));
        }

        min_val
    }

    fn bounds(&self) -> (nalgebra::Point3<T>, nalgebra::Point3<T>) {
        use cfd_mesh::application::delaunay::dim3::sdf::FiniteCylinderSdf;
        use nalgebra::{Point3, Vector3};
        
        let half = <T as cfd_mesh::domain::core::Scalar>::from_f64(2.0);
        
        let parent = FiniteCylinderSdf::new(
            Point3::new(-self.l_parent, T::zero(), T::zero()),
            Point3::new(T::zero(), T::zero(), T::zero()),
            self.d_parent / half,
        );
        let (mut min_pt, mut max_pt) = parent.bounds();
        
        for i in 0..3 {
            let angle = self.branching_angles[i];
            let dir = Vector3::new(num_traits::Float::cos(angle), num_traits::Float::sin(angle), T::zero());
            let daughter = FiniteCylinderSdf::new(
                Point3::new(T::zero(), T::zero(), T::zero()),
                Point3::new(T::zero(), T::zero(), T::zero()) + dir * self.l_daughters[i],
                self.d_daughters[i] / half,
            );
            let (d_min, d_max) = daughter.bounds();
            min_pt = Point3::new(
                num_traits::Float::min(min_pt.x, d_min.x),
                num_traits::Float::min(min_pt.y, d_min.y),
                num_traits::Float::min(min_pt.z, d_min.z),
            );
            max_pt = Point3::new(
                num_traits::Float::max(max_pt.x, d_max.x),
                num_traits::Float::max(max_pt.y, d_max.y),
                num_traits::Float::max(max_pt.z, d_max.z),
            );
        }
        
        (min_pt, max_pt)
    }
}
