use crate::topology::model::{BlueprintTopologySpec, SplitKind};

impl BlueprintTopologySpec {
    /// Human-readable topology name derived from the split-stage sequence.
    #[must_use]
    pub fn display_name(&self) -> String {
        if self.split_stages.is_empty() {
            return self.design_name.clone();
        }

        let stage_labels: Vec<String> = self
            .split_stages
            .iter()
            .map(|stage| match stage.split_kind {
                SplitKind::NFurcation(2) => "Bi".to_string(),
                SplitKind::NFurcation(3) => "Tri".to_string(),
                SplitKind::NFurcation(4) => "Quad".to_string(),
                SplitKind::NFurcation(5) => "Penta".to_string(),
                SplitKind::NFurcation(n) => format!("{n}Way"),
            })
            .collect();
        let tree_label = stage_labels.join("→");

        let suffix = if self.has_venturi() {
            format!(" + {}× Venturi", self.venturi_count())
        } else if self.has_serpentine() {
            let arm_count = self.serpentine_arm_count();
            format!(" + {arm_count}× Serpentine")
        } else {
            String::new()
        };

        format!("{tree_label}{suffix}")
    }

    /// Short (<=5 char) topology code for candidate IDs.
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

    /// Split-stage sequence label for render hints (e.g. `Tri→Bi`).
    #[must_use]
    pub fn stage_sequence_label(&self) -> String {
        if self.split_stages.is_empty() {
            return self.short_code();
        }
        self.split_stages
            .iter()
            .map(|stage| match stage.split_kind {
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

    /// Unique topology tag for candidate ID generation.
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
}
