//! Derived topology queries computed from [`BlueprintTopologySpec`] structure.
//!
//! These methods replace the per-variant `match` dispatch formerly in
//! `cfd-optim::DesignTopology`.  Every query is derived from the declarative
//! spec, not from an enum variant name — enabling the GA to compose arbitrary
//! topologies without extending an enum.

use super::model::{
    BlueprintTopologySpec, BranchRole, BranchSpec, SplitKind, TreatmentActuationMode,
};

// ── Venturi queries ──────────────────────────────────────────────────────────

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
            .map(|p| p.serial_throat_count as usize)
            .sum::<usize>()
            .max(if self.has_venturi() { 1 } else { 0 })
    }

    /// Number of venturi stages running in PARALLEL (for flow-per-venturi).
    ///
    /// For serial venturis on the same channel, all share the full flow,
    /// so they don't divide Q.  This returns the number of *distinct*
    /// target channels that carry venturi placements.
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
            .map(|p| p.target_channel_id.as_str())
            .collect();
        unique_channels.sort_unstable();
        unique_channels.dedup();
        unique_channels.len()
    }

    /// Maximum serial venturi stages on any single flow path.
    ///
    /// Scans all placements and returns the highest
    /// `serial_throat_count` value.  Returns 0 when no venturis exist.
    #[must_use]
    pub fn serial_venturi_stages(&self) -> usize {
        if self.has_series_path() && self.split_stages.is_empty() {
            return self
                .venturi_placements
                .iter()
                .map(|p| p.serial_throat_count as usize)
                .sum();
        }
        self.venturi_placements
            .iter()
            .map(|p| p.serial_throat_count as usize)
            .max()
            .unwrap_or(0)
    }

    /// Whether venturi treatment mode is active.
    #[must_use]
    pub fn uses_venturi_treatment(&self) -> bool {
        self.treatment_mode == TreatmentActuationMode::VenturiCavitation && self.has_venturi()
    }
}

// ── Serpentine / distribution queries ─────────────────────────────────────────

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
    ///
    /// Count of branches that carry a `SerpentineSpec` across all stages.
    /// For non-serpentine topologies returns 1 (the main flow path).
    #[must_use]
    pub fn serpentine_arm_count(&self) -> usize {
        let count: usize = self
            .split_stages
            .iter()
            .flat_map(|stage| stage.branches.iter())
            .map(|branch| branch.route.serpentine.is_some() as usize)
            .sum::<usize>()
            + self
                .series_channels
                .iter()
                .map(|channel| channel.route.serpentine.is_some() as usize)
                .sum::<usize>()
            + self
                .parallel_channels
                .iter()
                .map(|channel| channel.route.serpentine.is_some() as usize)
                .sum::<usize>();
        count.max(1)
    }

    /// Whether this topology includes a distribution network (any split stage).
    #[must_use]
    pub fn has_distribution(&self) -> bool {
        !self.split_stages.is_empty() || self.has_parallel_paths()
    }
}

// ── Branch / geometry queries ────────────────────────────────────────────────

impl BlueprintTopologySpec {
    /// Number of terminal (leaf) branches.
    ///
    /// Computed as the product of branch counts across all split stages.
    /// For a Bi→Tri tree this is 2 × 3 = 6 leaf channels.
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
    ///
    /// For selective routing trees, this is typically 1 per stage (the center arm).
    /// The total treatment leaf count for a multi-level tree where only the center
    /// arm continues is 1 (single treatment path through the entire tree),
    /// not the product of branch counts.
    #[must_use]
    pub fn treatment_branch_count(&self) -> usize {
        if self.split_stages.is_empty() {
            return self.treatment_channel_ids().len().max(1);
        }
        // Count the final-stage treatment branches only, since each stage
        // has exactly 1 treatment_path branch that feeds into the next stage.
        self.split_stages.last().map_or(1, |stage| {
            stage
                .branches
                .iter()
                .filter(|b| b.treatment_path)
                .count()
                .max(1)
        })
    }

    /// Number of physical outlet ports.
    ///
    /// All topologies produced by composite closed-loop presets converge to
    /// exactly one outlet.
    #[must_use]
    pub fn outlet_count(&self) -> usize {
        1
    }

    /// Fraction of the 36 treatment wells covered by this topology's channels.
    ///
    /// Derived from the topology structure:
    /// - Topologies with ≥ 2 split stages or serpentine arms → full coverage (1.0)
    /// - Single-point venturi with no distribution → 1/36
    /// - Cell separation (1 stage, 2+ branches, 1 treatment) → 0.5
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

        // Single-stage split: check if it's a cell-separation configuraton
        // (one treatment + bypass arms)
        let first = &self.split_stages[0];
        let n_treatment = first.branches.iter().filter(|b| b.treatment_path).count();
        let n_bypass = first.branches.len() - n_treatment;
        if n_bypass >= 1 && n_treatment == 1 {
            0.5 // center treatment + peripheral bypass
        } else {
            1.0 // fully distributed
        }
    }

    /// Human-readable topology name derived from the split-stage sequence.
    ///
    /// Examples: `"Tri→Bi + 2× Venturi"`, `"Serpentine Grid"`,
    /// `"Single Venturi"`.
    #[must_use]
    pub fn display_name(&self) -> String {
        if self.split_stages.is_empty() {
            return self.design_name.clone();
        }

        let stage_labels: Vec<String> = self
            .split_stages
            .iter()
            .map(|s| match s.split_kind {
                SplitKind::NFurcation(2) => "Bi".to_string(),
                SplitKind::NFurcation(3) => "Tri".to_string(),
                SplitKind::NFurcation(4) => "Quad".to_string(),
                SplitKind::NFurcation(5) => "Penta".to_string(),
                SplitKind::NFurcation(n) => format!("{n}Way"),
            })
            .collect();
        let tree_label = stage_labels.join("→");

        let leaf_count = self.terminal_branch_count();
        let suffix = if self.has_venturi() {
            format!(" + {leaf_count}× Venturi")
        } else if self.has_serpentine() {
            format!(" + {leaf_count}× Serpentine")
        } else {
            String::new()
        };

        format!("{tree_label}{suffix}")
    }

    /// Short (≤5 char) topology code for candidate IDs.
    ///
    /// Derived from split stages:
    /// - `"SV"` for single venturi
    /// - `"BV"` for Bi + venturi
    /// - `"TT"` for Tri→Tri
    /// - `"BTB"` for Bi→Tri→Bi
    #[must_use]
    pub fn short_code(&self) -> String {
        if self.split_stages.is_empty() {
            return if self.has_venturi() && self.has_serpentine() {
                "VSR".to_string()
            } else if self.has_venturi() && self.serial_venturi_stages() > 1 {
                "SDV".to_string()
            } else if self.has_venturi() {
                "SV".to_string()
            } else if self.has_serpentine() {
                "SERP".to_string()
            } else if self.has_parallel_paths() {
                "PARR".to_string()
            } else if self.has_series_path() {
                "LIN".to_string()
            } else {
                "SC".to_string()
            };
        }

        let mut code = String::with_capacity(8);
        for stage in &self.split_stages {
            code.push(match stage.split_kind {
                SplitKind::NFurcation(2) => 'B',
                SplitKind::NFurcation(3) => 'T',
                SplitKind::NFurcation(4) => 'Q',
                SplitKind::NFurcation(5) => 'P',
                SplitKind::NFurcation(_) => 'N',
            });
        }

        if self.has_venturi() {
            code.push('V');
        } else if self.has_serpentine() {
            code.push('S');
        }

        code
    }

    /// Split-stage sequence label for render hints (e.g. `"Tri→Bi"`).
    #[must_use]
    pub fn stage_sequence_label(&self) -> String {
        if self.split_stages.is_empty() {
            return self.short_code();
        }
        self.split_stages
            .iter()
            .map(|s| match s.split_kind {
                SplitKind::NFurcation(2) => "Bi".to_string(),
                SplitKind::NFurcation(3) => "Tri".to_string(),
                SplitKind::NFurcation(4) => "Quad".to_string(),
                SplitKind::NFurcation(5) => "Penta".to_string(),
                SplitKind::NFurcation(n) => format!("{n}Way"),
            })
            .collect::<Vec<_>>()
            .join("→")
    }

    /// Number of visible split layers for the topology legend.
    #[must_use]
    pub fn visible_split_layers(&self) -> usize {
        self.split_stages.len()
    }
}

// ── Classification queries (replace DesignTopology match arms) ───────────────

impl BlueprintTopologySpec {
    /// Number of split stages (depth of the split tree).
    ///
    /// Equivalent to the former `PrimitiveSplitSequence::levels()`.
    #[must_use]
    pub fn split_depth(&self) -> usize {
        self.split_stages.len()
    }

    /// Whether this topology is a leukapheresis-class design.
    ///
    /// Leukapheresis topologies operate at micro-scale channel dimensions
    /// (`D_h < 150 µm`, `κ_WBC > 0.15`) with diluted blood (`HCT ≈ 0.04`).
    /// Identified by having no venturi placements AND channel widths below
    /// the millifluidic threshold (or by checking the design_name convention).
    ///
    /// # Structural criterion
    ///
    /// A topology is leukapheresis-class when it has serpentine channels but
    /// no venturi placements and no treatment-path branches. These are
    /// the constriction-expansion, spiral, and parallel microchannel designs.
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
                .flat_map(|s| s.branches.iter())
                .any(|b| b.role == BranchRole::Treatment && b.treatment_path)
    }

    /// Whether the 2-population inertial separation model applies.
    ///
    /// True when the topology has a single-stage split with exactly one
    /// treatment branch and one or more bypass branches (the
    /// `CellSeparationVenturi` pattern). The Dean-flow margination focuses
    /// large stiff CTCs (17.5 µm, DI=0.15) to channel center while small
    /// deformable RBCs (7 µm, DI=0.85) migrate to the periphery.
    #[must_use]
    pub fn is_cell_separation(&self) -> bool {
        if self.split_stages.len() != 1 || !self.has_venturi() {
            return false;
        }
        let stage = &self.split_stages[0];
        let n_treat = stage.branches.iter().filter(|b| b.treatment_path).count();
        let n_bypass = stage
            .branches
            .iter()
            .filter(|b| b.role == BranchRole::RbcBypass)
            .count();
        stage.split_kind == SplitKind::NFurcation(2) && n_treat == 1 && n_bypass == 1
    }

    /// Whether the 3-population (cancer + WBC + RBC) separation model applies.
    ///
    /// True when any branch in the topology has a `WbcCollection` role,
    /// indicating 3-stream sorting: cancer treatment, WBC recovery, and
    /// RBC waste. Also true for topologies that have serpentine arms or
    /// complex trees with selective routing (trifurcation-first patterns).
    #[must_use]
    pub fn is_three_pop_separation(&self) -> bool {
        let has_wbc = self
            .split_stages
            .iter()
            .flat_map(|s| s.branches.iter())
            .any(|b| b.role == BranchRole::WbcCollection);
        let is_serpentine_based = self.has_serpentine() && !self.is_leukapheresis_topology();
        let is_selective_tree = self.is_selective_routing();
        has_wbc || is_serpentine_based || is_selective_tree
    }

    /// Whether this topology uses selective asymmetric routing.
    ///
    /// True when branch widths within any split stage are non-uniform
    /// (asymmetric splits for Zweifach-Fung cell sorting). Specifically,
    /// a stage with at least one `Treatment` branch and at least one
    /// non-treatment branch with a different width.
    #[must_use]
    pub fn is_selective_routing(&self) -> bool {
        for stage in &self.split_stages {
            let treatment_branches: Vec<&BranchSpec> =
                stage.branches.iter().filter(|b| b.treatment_path).collect();
            let bypass_branches: Vec<&BranchSpec> = stage
                .branches
                .iter()
                .filter(|b| !b.treatment_path)
                .collect();

            if !treatment_branches.is_empty() && !bypass_branches.is_empty() {
                let treat_w = treatment_branches[0].route.width_m;
                if bypass_branches
                    .iter()
                    .any(|b| (b.route.width_m - treat_w).abs() > 1e-12)
                {
                    return true;
                }
            }
        }
        false
    }

    /// Fraction of total flow that reaches the treatment zone.
    ///
    /// For symmetric splits all flow reaches treatment; for asymmetric
    /// selective routing, only the treatment-path fraction reaches the
    /// therapy channels.
    ///
    /// # Arguments
    ///
    /// * `center_frac_fn` — closure computing the hydraulic Q-fraction
    ///   from a geometric width fraction. For trifurcation with
    ///   Hagen-Poiseuille scaling: `tri_center_q_frac(width_frac)`.
    ///
    /// Falls back to width-ratio estimation when no closure is provided.
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
            let total_width: f64 = stage.branches.iter().map(|b| b.route.width_m).sum();
            if total_width <= 0.0 {
                continue;
            }
            let treatment_width: f64 = stage
                .branches
                .iter()
                .filter(|b| b.treatment_path)
                .map(|b| b.route.width_m)
                .sum();
            // For Hagen-Poiseuille in rectangular channels,
            // Q_branch / Q_total ≈ w_branch / w_total (first-order)
            frac *= treatment_width / total_width;
        }
        frac.clamp(0.0, 1.0)
    }

    /// Unique topology tag for candidate ID generation.
    ///
    /// Encodes the split sequence, venturi count, and serpentine status
    /// into a compact string suitable for embedding in candidate IDs.
    #[must_use]
    pub fn topology_tag(&self) -> String {
        let mut tag = self.short_code();
        let leaf_count = self.terminal_branch_count();
        if leaf_count > 1 {
            tag.push_str(&format!("-{leaf_count}"));
        }
        if self.serial_venturi_stages() > 1 {
            tag.push_str(&format!("-sv{}", self.serial_venturi_stages()));
        }
        tag
    }

    /// Whether the first split stage is a trifurcation.
    #[must_use]
    pub fn first_stage_is_trifurcation(&self) -> bool {
        self.split_stages
            .first()
            .map_or(false, |s| s.split_kind == SplitKind::NFurcation(3))
    }
}
