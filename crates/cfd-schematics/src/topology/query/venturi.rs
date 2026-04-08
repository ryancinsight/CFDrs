use crate::topology::model::{BlueprintTopologySpec, TreatmentActuationMode};

impl BlueprintTopologySpec {
    /// Whether this topology is expressed as a linear series channel path.
    #[must_use]
    pub fn has_series_path(&self) -> bool {
        !self.series_channels.is_empty()
    }

    /// Whether this topology is expressed as a parallel channel bundle.
    #[must_use]
    pub fn has_parallel_paths(&self) -> bool {
        !self.parallel_channels.is_empty()
    }

    /// Whether this topology includes at least one venturi placement.
    #[must_use]
    pub fn has_venturi(&self) -> bool {
        !self.venturi_placements.is_empty()
    }

    /// Total number of venturi throat placements (parallel + serial).
    #[must_use]
    pub fn venturi_count(&self) -> usize {
        self.venturi_placements
            .iter()
            .map(|placement| placement.serial_throat_count as usize)
            .sum::<usize>()
            .max(usize::from(self.has_venturi()))
    }

    /// Number of venturi stages running in parallel (for flow-per-venturi).
    ///
    /// For serial venturis on the same channel, all share the full flow,
    /// so they do not divide Q. This returns the number of distinct target
    /// channels that carry venturi placements.
    #[must_use]
    pub fn parallel_venturi_count(&self) -> usize {
        if self.venturi_placements.is_empty() {
            return 0;
        }
        if (self.has_series_path() || self.has_parallel_paths()) && self.split_stages.is_empty() {
            return 1;
        }
        let mut unique_channels: Vec<&str> = self
            .venturi_placements
            .iter()
            .map(|placement| placement.target_channel_id.as_str())
            .collect();
        unique_channels.sort_unstable();
        unique_channels.dedup();
        unique_channels.len()
    }

    /// Maximum serial venturi stages on any single flow path.
    #[must_use]
    pub fn serial_venturi_stages(&self) -> usize {
        if self.has_series_path() && self.split_stages.is_empty() {
            return self
                .venturi_placements
                .iter()
                .map(|placement| placement.serial_throat_count as usize)
                .sum();
        }
        self.venturi_placements
            .iter()
            .map(|placement| placement.serial_throat_count as usize)
            .max()
            .unwrap_or(0)
    }

    /// Whether venturi treatment mode is active.
    #[must_use]
    pub fn uses_venturi_treatment(&self) -> bool {
        self.treatment_mode == TreatmentActuationMode::VenturiCavitation && self.has_venturi()
    }
}