//! Primitive split-sequence catalog for selective-routing topologies.
//!
//! Each [`PrimitiveSplitSequence`] encodes a specific sequence of
//! bifurcation, trifurcation, quadfurcation, or pentafurcation stages used
//! to partition flow into treatment and bypass lanes.

use super::super::model::SplitKind;

/// Encodes one of the canonical primitive split-stage sequences used in the
/// Milestone 12 selective-routing design space.
///
/// # Nomenclature
///
/// Each variant name reads left-to-right as root → leaf.  `Bi`, `Tri`,
/// `Quad`, and `Penta` denote [`NFurcation(2)`], [`NFurcation(3)`],
/// [`NFurcation(4)`], and [`NFurcation(5)`] respectively.
///
/// # Example
///
/// ```
/// use cfd_schematics::topology::presets::PrimitiveSplitSequence;
///
/// let seq = PrimitiveSplitSequence::TriTri;
/// assert_eq!(seq.levels(), 2);
/// assert_eq!(seq.label(), "Tri→Tri");
/// assert_eq!(seq.to_split_kinds().len(), 2);
/// ```
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum PrimitiveSplitSequence {
    /// Single bifurcation (2 leaves).
    Bi,
    /// Single trifurcation (3 leaves).
    Tri,
    /// Trifurcation → bifurcation (6 leaves).
    TriBi,
    /// Trifurcation → trifurcation (9 leaves).
    TriTri,
    /// Trifurcation → bifurcation → bifurcation (12 leaves).
    TriBiBi,
    /// Trifurcation → bifurcation → trifurcation (18 leaves).
    TriBiTri,
    /// Trifurcation → trifurcation → bifurcation (18 leaves).
    TriTriBi,
    /// Trifurcation → trifurcation → trifurcation (27 leaves).
    TriTriTri,
    Quad,
    Penta,
    TriQuad,
    TriPenta,
    /// Quadfurcation → bifurcation (8 leaves).
    QuadBi,
    /// Quadfurcation → trifurcation (12 leaves).
    QuadTri,
    /// Quadfurcation → trifurcation → bifurcation (24 leaves).
    QuadTriBi,
    /// Pentafurcation → bifurcation (10 leaves).
    PentaBi,
    /// Pentafurcation → quadfurcation → bifurcation (40 leaves).
    PentaQuadBi,
    /// Pentafurcation → quadfurcation → trifurcation (60 leaves).
    PentaQuadTri,
    /// Pentafurcation → trifurcation (15 leaves).
    PentaTri,
    /// Pentafurcation → trifurcation → bifurcation (30 leaves).
    PentaTriBi,
}

impl PrimitiveSplitSequence {
    /// Number of split stages in this sequence.
    #[must_use]
    pub const fn levels(self) -> u8 {
        match self {
            Self::Bi | Self::Tri | Self::Quad | Self::Penta => 1,
            Self::TriBi | Self::TriTri | Self::TriQuad | Self::TriPenta
            | Self::QuadBi | Self::QuadTri | Self::PentaBi | Self::PentaTri => 2,
            Self::TriBiBi | Self::TriBiTri | Self::TriTriBi | Self::TriTriTri
            | Self::QuadTriBi | Self::PentaQuadBi | Self::PentaQuadTri | Self::PentaTriBi => 3,
        }
    }

    /// Human-readable label (e.g. `"Tri→Tri"`).
    #[must_use]
    pub const fn label(self) -> &'static str {
        match self {
            Self::Bi => "Bi",
            Self::Tri => "Tri",
            Self::TriBi => "Tri\u{2192}Bi",
            Self::TriTri => "Tri\u{2192}Tri",
            Self::TriBiBi => "Tri\u{2192}Bi\u{2192}Bi",
            Self::TriBiTri => "Tri\u{2192}Bi\u{2192}Tri",
            Self::TriTriBi => "Tri\u{2192}Tri\u{2192}Bi",
            Self::TriTriTri => "Tri\u{2192}Tri\u{2192}Tri",
            Self::Quad => "Quad",
            Self::Penta => "Penta",
            Self::TriQuad => "Tri\u{2192}Quad",
            Self::TriPenta => "Tri\u{2192}Penta",
            Self::QuadBi => "Quad\u{2192}Bi",
            Self::QuadTri => "Quad\u{2192}Tri",
            Self::QuadTriBi => "Quad\u{2192}Tri\u{2192}Bi",
            Self::PentaBi => "Penta\u{2192}Bi",
            Self::PentaQuadBi => "Penta\u{2192}Quad\u{2192}Bi",
            Self::PentaQuadTri => "Penta\u{2192}Quad\u{2192}Tri",
            Self::PentaTri => "Penta\u{2192}Tri",
            Self::PentaTriBi => "Penta\u{2192}Tri\u{2192}Bi",
        }
    }

    /// Convert this sequence to a [`Vec<SplitKind>`] for use with
    /// topology preset constructors.
    #[must_use]
    pub fn to_split_kinds(self) -> Vec<SplitKind> {
        self.label()
            .split('\u{2192}')
            .map(|segment| match segment {
                "Bi" => SplitKind::NFurcation(2),
                "Tri" => SplitKind::NFurcation(3),
                "Quad" => SplitKind::NFurcation(4),
                "Penta" => SplitKind::NFurcation(5),
                other => panic!("unknown primitive split segment '{other}'"),
            })
            .collect()
    }

    /// Returns `(has_intermediate_tri, has_any_tri, has_any_bi)`.
    ///
    /// Used by the parametric sweep to determine which fraction slices
    /// to iterate over when generating asymmetric-width candidates.
    #[must_use]
    pub fn metadata(self) -> (bool, bool, bool) {
        let label = self.label();
        let segments: Vec<&str> = label.split('\u{2192}').collect();
        let mut n_tri = 0;
        let mut has_bi = false;

        for segment in segments {
            if segment == "Tri" {
                n_tri += 1;
            }
            if segment == "Bi" {
                has_bi = true;
            }
        }

        (n_tri >= 2, n_tri > 0, has_bi)
    }
}

/// The tri-first primitive selective lineage sequences used by the
/// Milestone 12 design-space sweep.  Every entry starts with an
/// [`NFurcation(3)`] trifurcation to guarantee a center-lane treatment
/// path flanked by peripheral bypass channels.
pub const TRI_FIRST_SEQUENCES: [PrimitiveSplitSequence; 9] = [
    PrimitiveSplitSequence::Tri,
    PrimitiveSplitSequence::TriBi,
    PrimitiveSplitSequence::TriTri,
    PrimitiveSplitSequence::TriBiBi,
    PrimitiveSplitSequence::TriBiTri,
    PrimitiveSplitSequence::TriTriBi,
    PrimitiveSplitSequence::TriTriTri,
    PrimitiveSplitSequence::TriQuad,
    PrimitiveSplitSequence::TriPenta,
];

/// All primitive split sequences including non-tri-first variants for
/// broader parametric sweeps where the root stage may be Bi, Tri, Quad,
/// or Penta.
pub const ALL_SELECTIVE_SEQUENCES: [PrimitiveSplitSequence; 20] = [
    PrimitiveSplitSequence::Bi,
    PrimitiveSplitSequence::Tri,
    PrimitiveSplitSequence::TriBi,
    PrimitiveSplitSequence::TriTri,
    PrimitiveSplitSequence::TriBiBi,
    PrimitiveSplitSequence::TriBiTri,
    PrimitiveSplitSequence::TriTriBi,
    PrimitiveSplitSequence::TriTriTri,
    PrimitiveSplitSequence::Quad,
    PrimitiveSplitSequence::Penta,
    PrimitiveSplitSequence::TriQuad,
    PrimitiveSplitSequence::TriPenta,
    PrimitiveSplitSequence::QuadBi,
    PrimitiveSplitSequence::QuadTri,
    PrimitiveSplitSequence::QuadTriBi,
    PrimitiveSplitSequence::PentaBi,
    PrimitiveSplitSequence::PentaQuadBi,
    PrimitiveSplitSequence::PentaQuadTri,
    PrimitiveSplitSequence::PentaTri,
    PrimitiveSplitSequence::PentaTriBi,
];

/// The canonical Milestone 12 design-space sweep sequences.
///
/// Includes single-stage roots (Bi, Tri, Quad, Penta) plus multi-layer
/// Quad-first and Penta-first cascading split sequences.
pub const MILESTONE12_SWEEP_SEQUENCES: [PrimitiveSplitSequence; 12] = [
    PrimitiveSplitSequence::Bi,
    PrimitiveSplitSequence::Tri,
    PrimitiveSplitSequence::Quad,
    PrimitiveSplitSequence::Penta,
    PrimitiveSplitSequence::QuadBi,
    PrimitiveSplitSequence::QuadTri,
    PrimitiveSplitSequence::QuadTriBi,
    PrimitiveSplitSequence::PentaBi,
    PrimitiveSplitSequence::PentaQuadBi,
    PrimitiveSplitSequence::PentaQuadTri,
    PrimitiveSplitSequence::PentaTri,
    PrimitiveSplitSequence::PentaTriBi,
];

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn all_tri_first_sequences_start_with_trifurcation() {
        for &seq in &TRI_FIRST_SEQUENCES {
            let kinds = seq.to_split_kinds();
            assert_eq!(
                kinds[0],
                SplitKind::NFurcation(3),
                "{} must start with trifurcation",
                seq.label()
            );
        }
    }

    #[test]
    fn levels_match_split_kinds_length() {
        for &seq in &TRI_FIRST_SEQUENCES {
            assert_eq!(
                seq.levels() as usize,
                seq.to_split_kinds().len(),
                "{}: levels() disagrees with to_split_kinds().len()",
                seq.label()
            );
        }
    }

    #[test]
    fn metadata_consistency() {
        let (int_tri, any_tri, any_bi) = PrimitiveSplitSequence::Bi.metadata();
        assert!(!int_tri, "Bi has no tri stages");
        assert!(!any_tri, "Bi contains no tri");
        assert!(any_bi, "Bi contains bi");

        let (int_tri, any_tri, any_bi) = PrimitiveSplitSequence::TriTri.metadata();
        assert!(int_tri, "TriTri has ≥2 tri stages");
        assert!(any_tri, "TriTri contains tri");
        assert!(!any_bi, "TriTri has no bi stages");

        let (int_tri, any_tri, any_bi) = PrimitiveSplitSequence::TriBi.metadata();
        assert!(!int_tri, "TriBi has only 1 tri");
        assert!(any_tri, "TriBi contains tri");
        assert!(any_bi, "TriBi contains bi");
    }
}
