//! 3D Trifurcation geometry modeling
//!
//! Models a vessel branching into three daughter vessels.
//! Follows the same principles as Bifurcation3D but with three outlets.

use crate::bifurcation::ConicalTransition;
use cfd_core::conversion::SafeFromF64;
use nalgebra::{Point3, RealField, Vector3};
use num_traits::{FromPrimitive, ToPrimitive};
use serde::{Deserialize, Serialize};

/// 3D Trifurcation geometry
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TrifurcationGeometry3D<T: RealField + Copy> {
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

impl<T: RealField + Copy + FromPrimitive + ToPrimitive + SafeFromF64> TrifurcationGeometry3D<T> {
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
        let three = T::from_f64_or_one(3.0);
        let sum_daughters = self.d_daughters[0].powi(3)
            + self.d_daughters[1].powi(3)
            + self.d_daughters[2].powi(3);
        self.d_parent.powi(3) - sum_daughters
    }
}
