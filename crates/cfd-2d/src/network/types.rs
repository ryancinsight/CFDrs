use cfd_schematics::domain::model::CrossSectionSpec;
use cfd_schematics::domain::model::NetworkBlueprint;
use cfd_schematics::domain::therapy_metadata::TherapyZone;
use nalgebra::RealField;
use num_traits::{Float, FromPrimitive, ToPrimitive};

use super::projection::NetworkProjectionSummary;
use crate::solvers::ns_fvm::{NavierStokesSolver2D, SolveResult};

use super::{ChannelReferenceTrace, NetworkReferenceTrace};

/// Per-channel projection summary from schematics into the 2D solver domain.
#[derive(Debug, Clone)]
pub struct ChannelProjectionSummary<T> {
    /// Blueprint channel identifier.
    pub channel_id: String,
    /// Solver-domain length after projection [m].
    pub grid_length_m: T,
    /// Solver-domain width after projection [m].
    pub grid_width_m: T,
    /// Flattened path length used for projection [m].
    pub path_length_m: T,
    /// Schematic x-span of the routed path [m].
    pub path_span_x_m: T,
    /// Schematic y-span of the routed path [m].
    pub path_span_y_m: T,
    /// Number of fluid cells flagged by the projection.
    pub fluid_cell_count: usize,
    /// Fraction of solver cells occupied by fluid.
    pub fluid_fraction: T,
}

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
    /// Schematics-driven projection metadata for this channel.
    pub projection: ChannelProjectionSummary<T>,
    /// Estimated wall shear stress from the blueprint cross-section model [Pa].
    pub wall_shear_pa: T,
    /// Maximum wall shear stress extracted from the solved 2D field [Pa].
    pub field_wall_shear_max_pa: T,
    /// Mean wall shear stress extracted from the solved 2D field [Pa].
    pub field_wall_shear_mean_pa: T,
    /// Mean inlet pressure extracted from the solved 2D field [Pa].
    pub field_inlet_pressure_pa: T,
    /// Mean outlet pressure extracted from the solved 2D field [Pa].
    pub field_outlet_pressure_pa: T,
    /// Inlet-to-outlet pressure drop extracted from the solved 2D field [Pa].
    pub field_pressure_drop_pa: T,
    /// Effective hydraulic resistance extracted from the solved 2D field [Pa·s/m³].
    pub field_effective_resistance_pa_s_per_m3: T,
    /// Field-integrated outlet flow rate reconstructed from the solved 2D outlet profile [m^3/s].
    pub field_outlet_flow_m3_s: T,
    /// Difference between field-integrated outlet flow and the 1D reference flow [m^3/s].
    pub field_outlet_flow_error_m3_s: T,
    /// Relative outlet-flow error against the 1D reference flow [%].
    pub field_outlet_flow_error_pct: T,
    /// Transit time through the channel [s].
    pub transit_time_s: T,
    /// Eulerian-Lagrangian separation efficiency over the solved 2D field [%].
    pub field_separation_efficiency_pct: Option<T>,
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

/// Combined network result for the projected schematics-driven solve path.
#[derive(Debug, Clone)]
pub struct ProjectedNetwork2dResult<T> {
    /// Compatibility result returned by the existing per-channel solve path.
    pub result: Network2dResult<T>,
    /// Projection summaries for every channel in blueprint order.
    pub projection: NetworkProjectionSummary<T>,
}

/// Coupled network result returned by [`Network2DSolver::solve_coupled`].
#[derive(Debug, Clone)]
pub struct CoupledNetwork2dResult<T> {
    /// Final coupled 2D network result after the outer pressure-resistance iteration.
    pub result: Network2dResult<T>,
    /// Projection summaries for every channel in blueprint order.
    pub projection: NetworkProjectionSummary<T>,
    /// Number of outer coupled iterations executed.
    pub coupling_iterations: usize,
    /// Whether the coupled outer loop satisfied the requested tolerance.
    pub converged: bool,
    /// Maximum relative flow change observed in the final coupled iteration [%].
    pub max_relative_flow_change_pct: T,
    /// Maximum relative resistance change observed in the final coupled iteration [%].
    pub max_relative_resistance_change_pct: T,
}

/// Internal per-channel entry stored in [`Network2DSolver`].
pub(crate) struct Channel2dEntry<T: RealField + Copy + Float + FromPrimitive + ToPrimitive> {
    pub(crate) id: String,
    pub(crate) therapy_zone: TherapyZone,
    pub(crate) is_venturi_throat: bool,
    pub(crate) projection: ChannelProjectionSummary<T>,
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
    pub(crate) blueprint: NetworkBlueprint,
    pub(crate) channels: Vec<Channel2dEntry<T>>,
    pub(crate) reference_trace: NetworkReferenceTrace<T>,
    pub(crate) projection: NetworkProjectionSummary<T>,
    pub(crate) reference_density_kg_m3: f64,
    pub(crate) reference_viscosity_pa_s: f64,
    pub(crate) separation_tracking_enabled: bool,
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

    /// Return the cached per-channel projection summaries for the current solver state.
    #[must_use]
    pub fn projection_summary(&self) -> NetworkProjectionSummary<T> {
        self.projection.clone()
    }

    /// Borrow the cached projection summary without cloning.
    #[must_use]
    pub fn projection_summary_ref(&self) -> &NetworkProjectionSummary<T> {
        &self.projection
    }
}
