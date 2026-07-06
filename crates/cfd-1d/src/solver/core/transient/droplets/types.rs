use crate::scalar::Cfd1dScalar;
use crate::solver::core::transient::composition::MixtureComposition;
use eunomia::FloatElement;
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
pub struct DropletSplitPolicy<T: Cfd1dScalar + Copy> {
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

impl<T: Cfd1dScalar + Copy + FloatElement> Default for DropletSplitPolicy<T> {
    fn default() -> Self {
        Self {
            mode: SplitMode::AutoFlowWeighted,
            min_secondary_flow_fraction: <T as FloatElement>::from_f64(0.2),
            min_child_volume: <T as FloatElement>::from_f64(1e-15),
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
pub struct DropletInjection<T: Cfd1dScalar + Copy> {
    /// Unique droplet id.
    pub droplet_id: i32,
    /// Carrier fluid id (for provenance).
    pub fluid_id: i32,
    /// Droplet volume \[m³].
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
pub struct DropletBoundary<T: Cfd1dScalar + Copy> {
    /// Channel index.
    pub channel_index: usize,
    /// Relative position in channel [0, 1].
    pub relative_position: T,
}

/// Occupancy span for one channel segment [start, end].
#[derive(Debug, Clone)]
pub struct ChannelOccupancy<T: Cfd1dScalar + Copy> {
    /// Channel index.
    pub channel_index: usize,
    /// Start of occupied interval in channel coordinates.
    pub start: T,
    /// End of occupied interval in channel coordinates.
    pub end: T,
}

/// Position of a droplet while in network.
#[derive(Debug, Clone)]
pub struct DropletPosition<T: Cfd1dScalar + Copy> {
    /// Current edge index.
    pub channel_index: usize,
    /// Relative position in edge [0, 1].
    pub relative_position: T,
}

/// Per-timepoint droplet snapshot.
///
/// # Theorem - Finite-Length Occupancy Projection
///
/// For a droplet snapshot in [`DropletState::Network`], the occupied channel
/// set is the ordered unique projection of `occupancy_spans.channel_index`.
///
/// **Proof sketch**: The finite-length representation is the single
/// authoritative geometric state because each [`ChannelOccupancy`] stores an
/// occupied interval `[start, end]` in one channel. Projecting those intervals
/// onto channel indices and removing duplicates preserves every channel with
/// finite-length occupancy while discarding only interval extent. Since the
/// projection is computed on demand by [`DropletSnapshot::occupied_channels`],
/// there is no stored point-droplet channel set that can diverge from span
/// geometry.
#[derive(Debug, Clone)]
pub struct DropletSnapshot<T: Cfd1dScalar + Copy> {
    /// Droplet id.
    pub droplet_id: i32,
    /// Current state.
    pub state: DropletState,
    /// Current position when in network.
    pub position: Option<DropletPosition<T>>,
    /// Occupancy spans for finite-length tracking.
    pub occupancy_spans: Vec<ChannelOccupancy<T>>,
    /// Boundary points for finite-length tracking.
    pub boundaries: Vec<DropletBoundary<T>>,
    /// Total volume currently represented by this droplet \[m³].
    pub total_volume: T,
    /// Fluid id associated with droplet.
    pub fluid_id: i32,
    /// Optional local mixture sampled from composition pipeline.
    pub local_mixture: Option<MixtureComposition<T>>,
}

impl<T: Cfd1dScalar + Copy> DropletSnapshot<T> {
    /// Compute the ordered unique channel projection of finite-length spans.
    #[must_use]
    pub fn occupied_channels(&self) -> Vec<usize> {
        let mut channels = Vec::with_capacity(self.occupancy_spans.len());
        for span in &self.occupancy_spans {
            if !channels.contains(&span.channel_index) {
                channels.push(span.channel_index);
            }
        }
        channels
    }

    /// Backwards-compatible name for the finite-span projection.
    #[must_use]
    pub fn occupied_channels_from_spans(&self) -> Vec<usize> {
        self.occupied_channels()
    }

    /// Validate that finite-length occupancy spans define the snapshot.
    #[must_use]
    pub fn has_consistent_finite_length_occupancy(&self) -> bool {
        self.occupancy_spans
            .iter()
            .all(|span| span.start <= span.end && span.start >= T::zero() && span.end <= T::one())
    }
}

/// Droplet tracking state at one timepoint.
#[derive(Debug, Clone)]
pub struct DropletTrackingState<T: Cfd1dScalar + Copy> {
    /// Simulation time.
    pub time: T,
    /// Droplet snapshots keyed by droplet id.
    pub droplets: HashMap<i32, DropletSnapshot<T>>,
}

#[derive(Debug, Clone)]
pub(crate) struct ActiveDroplet<T: Cfd1dScalar + Copy> {
    pub(crate) state: DropletState,
    pub(crate) branches: Vec<DropletBranch<T>>,
}

#[derive(Debug, Clone)]
pub(crate) struct DropletBranch<T: Cfd1dScalar + Copy> {
    pub(crate) channel_index: usize,
    pub(crate) center: T,
    pub(crate) volume: T,
}

#[cfg(test)]
mod tests {
    use super::*;

    fn snapshot_with_spans(spans: Vec<ChannelOccupancy<f64>>) -> DropletSnapshot<f64> {
        DropletSnapshot {
            droplet_id: 7,
            state: DropletState::Network,
            position: None,
            occupancy_spans: spans,
            boundaries: Vec::new(),
            total_volume: 0.0,
            fluid_id: 1,
            local_mixture: None,
        }
    }

    #[test]
    fn occupied_channels_are_derived_from_finite_length_spans() {
        let snapshot = snapshot_with_spans(vec![
            ChannelOccupancy {
                channel_index: 3,
                start: 0.10,
                end: 0.30,
            },
            ChannelOccupancy {
                channel_index: 5,
                start: 0.00,
                end: 0.20,
            },
            ChannelOccupancy {
                channel_index: 3,
                start: 0.35,
                end: 0.60,
            },
        ]);

        assert_eq!(snapshot.occupied_channels(), vec![3, 5]);
        assert_eq!(
            snapshot.occupied_channels(),
            snapshot.occupied_channels_from_spans()
        );
        assert!(snapshot.has_consistent_finite_length_occupancy());
    }

    #[test]
    fn inconsistent_span_interval_is_rejected_by_snapshot_contract() {
        let snapshot = snapshot_with_spans(vec![ChannelOccupancy {
            channel_index: 2,
            start: 0.70,
            end: 0.20,
        }]);

        assert!(!snapshot.has_consistent_finite_length_occupancy());
    }
}
