use super::super::super::builder::primitive_selective_candidate;
use super::parameters::CandidateParams;
use crate::constraints::CHANNEL_HEIGHT_M;
use crate::domain::BlueprintCandidate;
use cfd_schematics::TreatmentActuationMode;

impl CandidateParams {
    /// Materialize the lightweight parameters into a full `BlueprintCandidate`
    /// containing the fully routed graph, polygons, and topology spec.
    #[must_use]
    pub fn materialize(&self) -> BlueprintCandidate {
        let seq_tag = self.seq_tag();
        let is_acoustic = self.treatment_actuation_mode == TreatmentActuationMode::UltrasoundOnly;

        let id = if is_acoustic {
            format!(
                "{:04}-PST-{}-pcf{}-tcf{}-btf{}-uo-q{:.0}ml-g{:.0}kPa-w{:.0}um-h{}-n{}",
                self.idx,
                seq_tag,
                (self.pretri_center_frac * 1000.0).round() as u32,
                (self.terminal_tri_center_frac * 1000.0).round() as u32,
                (self.bi_treat_frac * 1000.0).round() as u32,
                self.q * 6e7,
                self.gauge * 1e-3,
                self.w_ch * 1e6,
                (CHANNEL_HEIGHT_M * 1e6) as u32,
                self.n_segs,
            )
        } else {
            let tl_factor = if self.d_throat > 0.0 {
                self.throat_len / self.d_throat
            } else {
                0.0
            };
            format!(
                "{:04}-PST-{}-pcf{}-tcf{}-btf{}-vt{}-q{:.0}ml-g{:.0}kPa-dt{:.0}um-tl{}-w{:.0}um-h{}-n{}",
                self.idx,
                seq_tag,
                (self.pretri_center_frac * 1000.0).round() as u32,
                (self.terminal_tri_center_frac * 1000.0).round() as u32,
                (self.bi_treat_frac * 1000.0).round() as u32,
                self.vt_count,
                self.q * 6e7,
                self.gauge * 1e-3,
                self.d_throat * 1e6,
                tl_factor.round() as u32,
                self.w_ch * 1e6,
                (CHANNEL_HEIGHT_M * 1e6) as u32,
                self.n_segs,
            )
        };

        primitive_selective_candidate(
            id,
            &self.request,
            self.q,
            self.gauge,
            self.d_throat,
            self.throat_len,
            self.w_ch,
            self.n_segs,
            self.pretri_center_frac,
            self.terminal_tri_center_frac,
            self.bi_treat_frac,
            self.treatment_actuation_mode,
            self.vt_count,
        )
    }
}
