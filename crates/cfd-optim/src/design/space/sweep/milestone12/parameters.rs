use cfd_schematics::topology::presets::Milestone12TopologyRequest;
use cfd_schematics::TreatmentActuationMode;
use std::sync::Arc;

/// Lightweight parameter structure defining a Milestone 12 candidate without
/// allocating its multi-kilobyte blueprint geometry. ~100 bytes.
#[derive(Clone, Debug)]
pub struct CandidateParams {
    pub idx: u32,
    pub request: Arc<Milestone12TopologyRequest>,
    pub q: f64,
    pub gauge: f64,
    pub d_throat: f64,
    pub throat_len: f64,
    pub w_ch: f64,
    pub n_segs: usize,
    pub pretri_center_frac: f64,
    pub terminal_tri_center_frac: f64,
    pub bi_treat_frac: f64,
    pub treatment_actuation_mode: TreatmentActuationMode,
    pub vt_count: u8,
}

impl CandidateParams {
    /// Check if this candidate relies on venturi cavitation (Option 2).
    #[must_use]
    pub fn is_venturi(&self) -> bool {
        self.treatment_actuation_mode == TreatmentActuationMode::VenturiCavitation
            && self.vt_count > 0
    }

    /// Obtain the topological sequence tag (e.g. "Bi" or "TriTri").
    #[must_use]
    pub fn seq_tag(&self) -> String {
        self.request.design_name.replace('\u{2192}', "")
    }
}
