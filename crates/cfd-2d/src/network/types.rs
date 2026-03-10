use cfd_schematics::domain::model::CrossSectionSpec;
use cfd_schematics::domain::therapy_metadata::TherapyZone;
use nalgebra::RealField;
use num_traits::{Float, FromPrimitive, ToPrimitive};

use crate::solvers::ns_fvm::{NavierStokesSolver2D, SolveResult};

use super::{ChannelReferenceTrace, NetworkReferenceTrace};

/// Per-channel 2D solve result, including solver convergence data, haemolysis,
/// and outlet-flow agreement against the 1D reference channel solution.
#[derive(Debug, Clone)]
pub struct Channel2dResult<T> {
    /// Blueprint channel ID.
    pub channel_id: String,
    /// Therapy zone classification.
    pub therapy_zone: TherapyZone,
    /// Whether this channel contains a venturi throat geometry.
    pub is_venturi_throat: bool,
    /// SIMPLE solver convergence result.
    pub solve_result: SolveResult<T>,
    /// Estimated wall shear stress from the blueprint cross-section model [Pa].
    pub wall_shear_pa: T,
    /// Maximum wall shear stress extracted from the solved 2D field [Pa].
    pub field_wall_shear_max_pa: T,
    /// Mean wall shear stress extracted from the solved 2D field [Pa].
    pub field_wall_shear_mean_pa: T,
    /// Field-integrated outlet flow rate reconstructed from the solved 2D outlet profile [m^3/s].
    pub field_outlet_flow_m3_s: T,
    /// Difference between field-integrated outlet flow and the 1D reference flow [m^3/s].
    pub field_outlet_flow_error_m3_s: T,
    /// Relative outlet-flow error against the 1D reference flow [%].
    pub field_outlet_flow_error_pct: T,
    /// Transit time through the channel [s].
    pub transit_time_s: T,
    /// Giersiepen (1990) haemolysis index contribution (dimensionless).
    pub hemolysis_index: T,
    /// Authoritative cfd-1d reference trace for this channel.
    pub reference_trace: ChannelReferenceTrace<T>,
}

/// Results of a full-network 2D solve.
#[derive(Debug, Clone)]
pub struct Network2dResult<T> {
    /// Per-channel solve results in blueprint channel order.
    pub channels: Vec<Channel2dResult<T>>,
    /// Total haemolysis index across all channels (sum).
    pub total_hemolysis_index: T,
    /// Number of channels that converged.
    pub converged_count: usize,
    /// Maximum per-channel outlet-flow error against the 1D reference [%].
    pub max_field_outlet_flow_error_pct: T,
    /// Mean per-channel outlet-flow error against the 1D reference [%].
    pub mean_field_outlet_flow_error_pct: T,
    /// Authoritative cfd-1d reference solve used to configure the 2D run.
    pub reference_trace: NetworkReferenceTrace<T>,
}

/// Internal per-channel entry stored in [`Network2DSolver`].
pub(crate) struct Channel2dEntry<T: RealField + Copy + Float + FromPrimitive + ToPrimitive> {
    pub(crate) id: String,
    pub(crate) therapy_zone: TherapyZone,
    pub(crate) is_venturi_throat: bool,
    /// Flow rate through this channel [m³/s].
    pub(crate) flow_rate_m3_s: f64,
    pub(crate) cross_section: CrossSectionSpec,
    pub(crate) cross_section_area_m2: f64,
    pub(crate) length_m: f64,
    pub(crate) viscosity_pa_s: f64,
    pub(crate) reference_trace: ChannelReferenceTrace<T>,
    pub(crate) solver: NavierStokesSolver2D<T>,
}

/// Full-network 2D solver: one configured 2D channel domain per blueprint channel.
pub struct Network2DSolver<T: RealField + Copy + Float + FromPrimitive + ToPrimitive> {
    pub(crate) channels: Vec<Channel2dEntry<T>>,
    pub(crate) reference_trace: NetworkReferenceTrace<T>,
}

impl<T> Network2DSolver<T>
where
    T: RealField + Copy + Float + FromPrimitive + ToPrimitive + std::fmt::Debug,
{
    /// Return the authoritative cfd-1d reference trace used to configure the 2D network.
    #[must_use]
    pub fn reference_trace(&self) -> &NetworkReferenceTrace<T> {
        &self.reference_trace
    }
}
