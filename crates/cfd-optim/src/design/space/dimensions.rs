//! Milestone 12 selective sweep dimensions.

use crate::constraints::TRIFURCATION_CENTER_FRACS;

pub(super) struct CavitationDimensions {
    pub flows: &'static [f64],
    pub gauges: &'static [f64],
    pub throats: &'static [f64],
}

static CAV_FLOWS: [f64; 3] = [5.000e-6, 6.667e-6, 8.333e-6];
static CAV_GAUGES: [f64; 3] = [300_000.0, 400_000.0, 500_000.0];
static CAV_THROATS: [f64; 3] = [30e-6, 35e-6, 45e-6];

impl CavitationDimensions {
    /// Return the active cavitation-tier dimensions.
    #[must_use]
    pub fn active() -> Self {
        #[cfg(any(test, debug_assertions))]
        {
            Self {
                flows: &CAV_FLOWS[0..1],
                gauges: &CAV_GAUGES[0..1],
                throats: &CAV_THROATS[0..1],
            }
        }
        #[cfg(not(any(test, debug_assertions)))]
        {
            Self {
                flows: &CAV_FLOWS,
                gauges: &CAV_GAUGES,
                throats: &CAV_THROATS,
            }
        }
    }
}

#[cfg(any(test, debug_assertions))]
fn primary_slice<T>(values: &'static [T]) -> &'static [T] {
    &values[..1]
}

#[cfg(not(any(test, debug_assertions)))]
fn primary_slice<T>(values: &'static [T]) -> &'static [T] {
    values
}

pub(super) const SELECTIVE_TRI_CENTER_FRACS: [f64; 3] = [0.33, 0.45, 0.55];
pub(super) const SELECTIVE_BI_TREAT_FRACS: [f64; 3] = [0.60, 0.68, 0.76];
pub(super) const SELECTIVE_CENTERLINE_VT_COUNTS: [u8; 2] = [1, 2];

static MILESTONE12_FLOWS: [f64; 4] = [1.333e-6, 1.667e-6, 2.000e-6, 3.333e-6];
static MILESTONE12_GAUGES: [f64; 4] = [25_000.0, 50_000.0, 100_000.0, 200_000.0];
static MILESTONE12_THROATS: [f64; 6] = [35e-6, 45e-6, 55e-6, 75e-6, 100e-6, 120e-6];
static MILESTONE12_WIDTHS: [f64; 3] = [4.0e-3, 6.0e-3, 8.0e-3];

pub(super) struct Milestone12Dimensions {
    pub flows: &'static [f64],
    pub gauges: &'static [f64],
    pub throats: &'static [f64],
    pub widths: &'static [f64],
    pub venturi_counts: &'static [u8],
    pub cavitation: CavitationDimensions,
}

/// Compute the appropriate fraction slices for a PST sequence's split-type
/// metadata.
///
/// Returns `(pretri_fracs, tri_center_fracs, bi_treat_fracs)` slices whose
/// lengths depend on whether the sequence contains intermediate tri-splits,
/// any tri-splits, or any bi-splits.
pub(super) fn pst_frac_slices(
    has_intermediate_tri: bool,
    has_any_tri: bool,
    has_any_bi: bool,
) -> (&'static [f64], &'static [f64], &'static [f64]) {
    static DEFAULT_PRETRI: [f64; 1] = [0.33];
    static DEFAULT_TRI: [f64; 1] = [1.0 / 3.0];
    static DEFAULT_BI: [f64; 1] = [0.68];

    let pretri = if has_intermediate_tri {
        &SELECTIVE_TRI_CENTER_FRACS[..]
    } else {
        &DEFAULT_PRETRI[..]
    };
    let tri = if has_any_tri {
        &TRIFURCATION_CENTER_FRACS
    } else {
        &DEFAULT_TRI[..]
    };
    let bi = if has_any_bi {
        &SELECTIVE_BI_TREAT_FRACS[..]
    } else {
        &DEFAULT_BI[..]
    };
    (pretri, tri, bi)
}

impl Milestone12Dimensions {
    #[must_use]
    pub fn active() -> Self {
        Self {
            flows: primary_slice(&MILESTONE12_FLOWS),
            gauges: primary_slice(&MILESTONE12_GAUGES),
            throats: primary_slice(&MILESTONE12_THROATS),
            widths: primary_slice(&MILESTONE12_WIDTHS),
            venturi_counts: primary_slice(&SELECTIVE_CENTERLINE_VT_COUNTS),
            cavitation: CavitationDimensions::active(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::Milestone12Dimensions;

    #[test]
    fn milestone12_dimensions_use_reduced_test_grid() {
        let dims = Milestone12Dimensions::active();
        assert_eq!(dims.flows.len(), 1);
        assert_eq!(dims.gauges.len(), 1);
        assert_eq!(dims.throats.len(), 1);
        assert_eq!(dims.widths.len(), 1);
        assert_eq!(dims.venturi_counts.len(), 1);
        assert_eq!(dims.cavitation.flows.len(), 1);
        assert_eq!(dims.cavitation.gauges.len(), 1);
        assert_eq!(dims.cavitation.throats.len(), 1);
    }
}
