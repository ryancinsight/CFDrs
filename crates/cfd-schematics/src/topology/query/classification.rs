use crate::topology::model::{BlueprintTopologySpec, BranchRole, BranchSpec, SplitKind};

impl BlueprintTopologySpec {
    /// Number of split stages (depth of the split tree).
    #[must_use]
    pub fn split_depth(&self) -> usize {
        self.split_stages.len()
    }

    /// Whether this topology is a leukapheresis-class design.
    #[must_use]
    pub fn is_leukapheresis_topology(&self) -> bool {
        if self.has_parallel_paths() && self.split_stages.is_empty() && !self.has_venturi() {
            return true;
        }
        if self.has_series_path() && self.split_stages.is_empty() && !self.has_venturi() {
            return self.series_channels.iter().any(|channel| {
                channel.channel_id.starts_with("narrow_")
                    || channel.channel_id.starts_with("spiral_")
                    || channel.channel_id.starts_with("ch_")
            });
        }
        !self.has_venturi()
            && self.has_serpentine()
            && !self
                .split_stages
                .iter()
                .flat_map(|stage| stage.branches.iter())
                .any(|branch| branch.role == BranchRole::Treatment && branch.treatment_path)
    }

    /// Whether the 2-population inertial separation model applies.
    #[must_use]
    pub fn is_cell_separation(&self) -> bool {
        if self.split_stages.len() != 1 || !self.has_venturi() {
            return false;
        }
        let stage = &self.split_stages[0];
        let n_treat = stage
            .branches
            .iter()
            .filter(|branch| branch.treatment_path)
            .count();
        let n_bypass = stage
            .branches
            .iter()
            .filter(|branch| branch.role == BranchRole::RbcBypass)
            .count();
        stage.split_kind == SplitKind::NFurcation(2) && n_treat == 1 && n_bypass == 1
    }

    /// Whether the 3-population (cancer + WBC + RBC) separation model applies.
    #[must_use]
    pub fn is_three_pop_separation(&self) -> bool {
        let has_wbc = self
            .split_stages
            .iter()
            .flat_map(|stage| stage.branches.iter())
            .any(|branch| branch.role == BranchRole::WbcCollection);
        let is_serpentine_based = self.has_serpentine() && !self.is_leukapheresis_topology();
        let is_selective_tree = self.is_selective_routing();
        has_wbc || is_serpentine_based || is_selective_tree
    }

    /// Whether this topology uses selective asymmetric routing.
    #[must_use]
    pub fn is_selective_routing(&self) -> bool {
        for stage in &self.split_stages {
            let treatment_branches: Vec<&BranchSpec> = stage
                .branches
                .iter()
                .filter(|branch| branch.treatment_path)
                .collect();
            let bypass_branches: Vec<&BranchSpec> = stage
                .branches
                .iter()
                .filter(|branch| !branch.treatment_path)
                .collect();

            if !treatment_branches.is_empty() && !bypass_branches.is_empty() {
                let treat_w = treatment_branches[0].route.width_m;
                if bypass_branches
                    .iter()
                    .any(|branch| (branch.route.width_m - treat_w).abs() > 1e-12)
                {
                    return true;
                }
            }
        }
        false
    }

    /// Fraction of total flow that reaches the treatment zone.
    #[must_use]
    pub fn therapy_channel_fraction(&self) -> f64 {
        if self.split_stages.is_empty() {
            return if self.treatment_channel_ids().is_empty() {
                0.0
            } else {
                1.0
            };
        }

        let mut frac = 1.0_f64;
        for stage in &self.split_stages {
            let total_width: f64 = stage
                .branches
                .iter()
                .map(|branch| branch.route.width_m)
                .sum();
            if total_width <= 0.0 {
                continue;
            }
            let treatment_width: f64 = stage
                .branches
                .iter()
                .filter(|branch| branch.treatment_path)
                .map(|branch| branch.route.width_m)
                .sum();
            frac *= treatment_width / total_width;
        }
        frac.clamp(0.0, 1.0)
    }

    /// Whether the first split stage is a trifurcation.
    #[must_use]
    pub fn first_stage_is_trifurcation(&self) -> bool {
        self.split_stages
            .first()
            .is_some_and(|stage| stage.split_kind == SplitKind::NFurcation(3))
    }
}
