use super::analytical::BernoulliVenturi;
use crate::scalar::Cfd2dScalar;
use crate::scalar::{self, from_f64};
use eunomia::{FloatElement, NumericElement};
use serde::{Deserialize, Serialize};

/// Solution to the Venturi flow problem
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct VenturiFlowSolution<T: Cfd2dScalar + Copy> {
    /// Inlet velocity \[m/s]
    pub u_inlet: T,
    /// Inlet pressure \[Pa]
    pub p_inlet: T,
    /// Maximum throat velocity \[m/s] (centerline of parabolic profile)
    pub u_throat: T,
    /// Area-averaged throat velocity \[m/s].
    pub u_throat_mean: T,
    /// Throat pressure \[Pa]
    pub p_throat: T,
    /// Outlet velocity \[m/s]
    pub u_outlet: T,
    /// Outlet pressure \[Pa]
    pub p_outlet: T,
    /// Pressure drop in throat \[Pa]
    pub dp_throat: T,
    /// Pressure recovery (outlet - inlet) \[Pa]
    pub dp_recovery: T,
    /// Pressure coefficient at throat
    pub cp_throat: T,
    /// Pressure recovery coefficient at outlet
    pub cp_recovery: T,
    /// Whether the iterative solver converged (always true for Bernoulli)
    pub converged: bool,
    /// SIMPLE iterations consumed (0 for the analytical Bernoulli path).
    pub iterations: usize,
    /// Final continuity residual reported by the iterative solver.
    pub final_residual: T,
}

impl<T: Cfd2dScalar + Copy + FloatElement> VenturiFlowSolution<T> {
    /// Create Venturi solution from Bernoulli
    pub fn from_bernoulli(bernoulli: &BernoulliVenturi<T>, p_outlet: T) -> Self {
        let u_throat = bernoulli.velocity_throat();
        let p_throat = bernoulli.pressure_throat();
        let cp_throat = bernoulli.pressure_coefficient_throat();

        let one_half = from_f64::<T>(0.5);
        let dynamic_pressure = one_half * bernoulli.rho * bernoulli.u_inlet * bernoulli.u_inlet;
        let cp_recovery = (p_outlet - bernoulli.p_inlet)
            / <T as NumericElement>::max_scalar(dynamic_pressure, from_f64::<T>(1.0));

        Self {
            u_inlet: bernoulli.u_inlet,
            p_inlet: bernoulli.p_inlet,
            u_throat,
            u_throat_mean: u_throat,
            p_throat,
            u_outlet: bernoulli.u_inlet,
            p_outlet,
            dp_throat: p_throat - bernoulli.p_inlet,
            dp_recovery: p_outlet - bernoulli.p_inlet,
            cp_throat,
            cp_recovery,
            converged: true,
            iterations: 0,
            final_residual: scalar::zero::<T>(),
        }
    }

    /// Verify energy conservation (should account for dissipation)
    pub fn energy_dissipation(&self, rho: T) -> T {
        let one_half = from_f64::<T>(0.5);

        let energy_inlet = self.p_inlet + one_half * rho * self.u_inlet * self.u_inlet;
        let energy_outlet = self.p_outlet + one_half * rho * self.u_outlet * self.u_outlet;

        <T as NumericElement>::max_scalar(energy_inlet - energy_outlet, scalar::zero::<T>())
    }
}
