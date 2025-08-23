//! Flow regime classification
//!
//! Provides flow regime identification based on dimensionless numbers
//! and flow characteristics.

use nalgebra::RealField;
use num_traits::ToPrimitive;
use serde::{Deserialize, Serialize};

/// Flow regime classification
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub enum FlowRegime {
    /// Stokes flow (Re << 1)
    Stokes,
    /// Laminar flow
    Laminar,
    /// Transitional flow
    Transitional,
    /// Turbulent flow
    Turbulent,
    /// Hypersonic flow
    Hypersonic,
}

/// Flow classifier based on dimensionless numbers
pub struct FlowClassifier;

impl FlowClassifier {
    /// Classify flow regime based on Reynolds number
    pub fn classify_by_reynolds<T: RealField + Copy + ToPrimitive>(reynolds: T) -> FlowRegime {
        let re = reynolds.to_f64().unwrap_or(0.0);
        
        if re < 0.1 {
            FlowRegime::Stokes
        } else if re < crate::constants::physics::dimensionless::reynolds::PIPE_CRITICAL_LOWER {
            FlowRegime::Laminar
        } else if re < crate::constants::physics::dimensionless::reynolds::PIPE_CRITICAL_UPPER {
            FlowRegime::Transitional
        } else {
            FlowRegime::Turbulent
        }
    }
    
    /// Classify flow regime based on Mach number
    pub fn classify_by_mach<T: RealField + Copy + ToPrimitive>(mach: T) -> FlowRegime {
        let ma = mach.to_f64().unwrap_or(0.0);
        
        if ma < crate::constants::physics::dimensionless::mach::INCOMPRESSIBLE_LIMIT {
            // Incompressible flow - further classify by Reynolds
            FlowRegime::Laminar // Would need Reynolds for proper classification
        } else if ma < crate::constants::physics::dimensionless::mach::TRANSONIC_LOWER {
            // Subsonic compressible
            FlowRegime::Turbulent // Typically turbulent at high speeds
        } else if ma < crate::constants::physics::dimensionless::mach::HYPERSONIC {
            // Supersonic
            FlowRegime::Turbulent
        } else {
            // Hypersonic
            FlowRegime::Hypersonic
        }
    }
    
    /// Classify based on multiple dimensionless numbers
    pub fn classify<T: RealField + Copy + ToPrimitive>(
        reynolds: Option<T>,
        mach: Option<T>,
        froude: Option<T>,
    ) -> FlowRegime {
        // Priority: Mach number for compressibility, then Reynolds for turbulence
        if let Some(ma) = mach {
            if ma.to_f64().unwrap_or(0.0) >= crate::constants::physics::dimensionless::mach::HYPERSONIC {
                return FlowRegime::Hypersonic;
            }
        }
        
        if let Some(re) = reynolds {
            return Self::classify_by_reynolds(re);
        }
        
        // Default to laminar if no information available
        FlowRegime::Laminar
    }
}

impl FlowRegime {
    /// Check if flow is viscous-dominated
    #[must_use] pub fn is_viscous_dominated(&self) -> bool {
        matches!(self, FlowRegime::Stokes | FlowRegime::Laminar)
    }
    
    /// Check if flow is inertia-dominated
    #[must_use] pub fn is_inertia_dominated(&self) -> bool {
        matches!(self, FlowRegime::Turbulent | FlowRegime::Hypersonic)
    }
    
    /// Check if flow requires turbulence modeling
    #[must_use] pub fn requires_turbulence_model(&self) -> bool {
        matches!(self, FlowRegime::Transitional | FlowRegime::Turbulent)
    }
    
    /// Get typical CFL number for this regime
    #[must_use] pub fn typical_cfl(&self) -> f64 {
        match self {
            FlowRegime::Stokes => 10.0,  // Can use large time steps
            FlowRegime::Laminar => 1.0,
            FlowRegime::Transitional => 0.5,
            FlowRegime::Turbulent => 0.3,
            FlowRegime::Hypersonic => 0.1,  // Need small time steps
        }
    }
}