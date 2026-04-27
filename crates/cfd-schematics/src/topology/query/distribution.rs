use crate::topology::model::BlueprintTopologySpec;

impl BlueprintTopologySpec {
    /// Whether any branch in this topology includes a serpentine path.
    #[must_use]
    pub fn has_serpentine(&self) -> bool {
        self.split_stages.iter().any(|stage| {
            stage
                .branches
                .iter()
                .any(|branch| branch.route.serpentine.is_some())
        }) || self
            .series_channels
            .iter()
            .any(|channel| channel.route.serpentine.is_some())
            || self
                .parallel_channels
                .iter()
                .any(|channel| channel.route.serpentine.is_some())
    }

    /// Number of parallel serpentine arms.
    #[must_use]
    pub fn serpentine_arm_count(&self) -> usize {
        let count: usize = self
            .split_stages
            .iter()
            .flat_map(|stage| stage.branches.iter())
            .map(|branch| usize::from(branch.route.serpentine.is_some()))
            .sum::<usize>()
            + self
                .series_channels
                .iter()
                .map(|channel| usize::from(channel.route.serpentine.is_some()))
                .sum::<usize>()
            + self
                .parallel_channels
                .iter()
                .map(|channel| usize::from(channel.route.serpentine.is_some()))
                .sum::<usize>();
        count.max(1)
    }

    /// Whether this topology includes a distribution network (any split stage).
    #[must_use]
    pub fn has_distribution(&self) -> bool {
        !self.split_stages.is_empty() || self.has_parallel_paths()
    }

    /// Number of terminal (leaf) branches.
    #[must_use]
    pub fn terminal_branch_count(&self) -> usize {
        if self.split_stages.is_empty() {
            if self.has_parallel_paths() {
                return self.parallel_channels.len().max(1);
            }
            return 1;
        }
        self.split_stages
            .iter()
            .map(|stage| stage.split_kind.branch_count())
            .product()
    }

    /// Number of treatment branches (branches with `treatment_path = true`).
    #[must_use]
    pub fn treatment_branch_count(&self) -> usize {
        if self.split_stages.is_empty() {
            return self.treatment_channel_ids().len().max(1);
        }
        self.split_stages.last().map_or(1, |stage| {
            stage
                .branches
                .iter()
                .filter(|branch| branch.treatment_path)
                .count()
                .max(1)
        })
    }

    /// Number of physical outlet ports.
    #[must_use]
    pub fn outlet_count(&self) -> usize {
        1
    }

    /// Fraction of the 36 treatment wells covered by this topology's channels.
    #[must_use]
    pub fn nominal_well_coverage(&self) -> f64 {
        const TREATMENT_WELL_COUNT: f64 = 36.0;

        if self.split_stages.is_empty() {
            if self.has_serpentine() {
                return 1.0;
            }
            let treatment_count = self.treatment_channel_ids().len();
            if self.has_parallel_paths() && treatment_count > 0 {
                return 1.0;
            }
            return if self.has_venturi() && treatment_count > 0 {
                (treatment_count as f64 / TREATMENT_WELL_COUNT)
                    .clamp(1.0 / TREATMENT_WELL_COUNT, 1.0)
            } else {
                1.0
            };
        }

        if self.has_serpentine() || self.split_stages.len() >= 2 {
            return 1.0;
        }

        let first = &self.split_stages[0];
        let n_treatment = first
            .branches
            .iter()
            .filter(|branch| branch.treatment_path)
            .count();
        let n_bypass = first.branches.len() - n_treatment;
        if n_bypass >= 1 && n_treatment == 1 {
            0.5
        } else {
            1.0
        }
    }
}
