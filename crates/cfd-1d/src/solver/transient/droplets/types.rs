use crate::solver::transient::composition::MixtureComposition;
use nalgebra::RealField;
use num_traits::FromPrimitive;
use std::collections::HashMap;

/// Split mode selection at junctions.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SplitMode {
    /// Use automatic flow/volume-aware split decision.
    AutoFlowWeighted,
    /// Force splitting across all outgoing branches.
    AlwaysSplit,
    /// Disable splitting and route to dominant branch only.
    NeverSplit,
}

/// Policy controlling droplet split behavior at junctions.
#[derive(Debug, Clone)]
pub struct DropletSplitPolicy<T: RealField + Copy> {
    /// Split decision mode.
    pub mode: SplitMode,
    /// Minimum secondary-flow fraction needed to trigger split in auto mode.
    /// Example: 0.2 means the non-dominant outgoing flow must be at least 20%.
    pub min_secondary_flow_fraction: T,
    /// Minimum child droplet volume allowed after splitting.
    pub min_child_volume: T,
    /// Maximum number of outgoing branches to distribute into.
    pub max_split_branches: usize,
}

impl<T: RealField + Copy + FromPrimitive> Default for DropletSplitPolicy<T> {
    fn default() -> Self {
        Self {
            mode: SplitMode::AutoFlowWeighted,
            min_secondary_flow_fraction: T::from_f64(0.2).unwrap_or_else(T::zero),
            min_child_volume: T::from_f64(1e-15).unwrap_or_else(T::zero),
            max_split_branches: 2,
        }
    }
}

/// High-level droplet state in the network lifecycle.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DropletState {
    /// Droplet is defined but not yet injected.
    Injection,
    /// Droplet currently travels in the network.
    Network,
    /// Droplet reached an outlet/sink.
    Sink,
    /// Droplet cannot progress due to missing valid outgoing path.
    Trapped,
}

/// Droplet injection definition.
#[derive(Debug, Clone)]
pub struct DropletInjection<T: RealField + Copy> {
    /// Unique droplet id.
    pub droplet_id: i32,
    /// Carrier fluid id (for provenance).
    pub fluid_id: i32,
    /// Droplet volume [m³].
    pub volume: T,
    /// Injection time.
    pub injection_time: T,
    /// Injection edge index.
    pub channel_index: usize,
    /// Relative position in edge [0, 1].
    pub relative_position: T,
}

/// Droplet boundary point on a channel.
#[derive(Debug, Clone)]
pub struct DropletBoundary<T: RealField + Copy> {
    /// Channel index.
    pub channel_index: usize,
    /// Relative position in channel [0, 1].
    pub relative_position: T,
}

/// Occupancy span for one channel segment [start, end].
#[derive(Debug, Clone)]
pub struct ChannelOccupancy<T: RealField + Copy> {
    /// Channel index.
    pub channel_index: usize,
    /// Start of occupied interval in channel coordinates.
    pub start: T,
    /// End of occupied interval in channel coordinates.
    pub end: T,
}

/// Position of a droplet while in network.
#[derive(Debug, Clone)]
pub struct DropletPosition<T: RealField + Copy> {
    /// Current edge index.
    pub channel_index: usize,
    /// Relative position in edge [0, 1].
    pub relative_position: T,
}

/// Per-timepoint droplet snapshot.
#[derive(Debug, Clone)]
pub struct DropletSnapshot<T: RealField + Copy> {
    /// Droplet id.
    pub droplet_id: i32,
    /// Current state.
    pub state: DropletState,
    /// Current position when in network.
    pub position: Option<DropletPosition<T>>,
    /// Occupied channel ids (point-droplet approximation uses 0 or 1 edge).
    pub occupied_channels: Vec<usize>,
    /// Occupancy spans for finite-length tracking.
    pub occupancy_spans: Vec<ChannelOccupancy<T>>,
    /// Boundary points for finite-length tracking.
    pub boundaries: Vec<DropletBoundary<T>>,
    /// Total volume currently represented by this droplet [m³].
    pub total_volume: T,
    /// Fluid id associated with droplet.
    pub fluid_id: i32,
    /// Optional local mixture sampled from composition pipeline.
    pub local_mixture: Option<MixtureComposition<T>>,
}

/// Droplet tracking state at one timepoint.
#[derive(Debug, Clone)]
pub struct DropletTrackingState<T: RealField + Copy> {
    /// Simulation time.
    pub time: T,
    /// Droplet snapshots keyed by droplet id.
    pub droplets: HashMap<i32, DropletSnapshot<T>>,
}

#[derive(Debug, Clone)]
pub(crate) struct ActiveDroplet<T: RealField + Copy> {
    pub(crate) state: DropletState,
    pub(crate) branches: Vec<DropletBranch<T>>,
}

#[derive(Debug, Clone)]
pub(crate) struct DropletBranch<T: RealField + Copy> {
    pub(crate) channel_index: usize,
    pub(crate) center: T,
    pub(crate) volume: T,
}
