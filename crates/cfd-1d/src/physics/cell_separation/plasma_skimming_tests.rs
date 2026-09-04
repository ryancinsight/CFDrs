//! Tests for the Pries phase-separation and plasma-skimming correlations.
//!
//! A sidecar in the form this workspace already uses (`chebyshev_tests.rs`):
//! the module was two fifths of the file it was written in.

use super::*;

const HT_NORMAL: f64 = 0.45;
const D_FEED: Length = Length::from_base(100.0e-6);

fn diameter_um(value: f64) -> Length {
    Length::from_base(value * 1.0e-6)
}

/// Equal 50/50 flow split with equal diameters: daughter hematocrit
/// should be approximately equal to feed hematocrit (symmetric case,
/// no skimming bias).
#[test]
fn test_plasma_skimming_equal_split() -> cfd_core::error::Result<()> {
    let ht = plasma_skimming_hematocrit(HT_NORMAL, 0.5, D_FEED, D_FEED)?;
    let rel_err = (ht - HT_NORMAL).abs() / HT_NORMAL;
    assert!(
        rel_err < 0.05,
        "Equal split hematocrit ({ht:.4}) should be ~feed ({HT_NORMAL:.4}), rel_err = {rel_err:.4}"
    );
    Ok(())
}

/// The smaller daughter branch (receiving 20% of flow) should get
/// less hematocrit than the feed due to plasma skimming.
#[test]
fn test_plasma_skimming_smaller_daughter_less_hematocrit() -> cfd_core::error::Result<()> {
    let ht_small = plasma_skimming_hematocrit(HT_NORMAL, 0.2, diameter_um(60.0), D_FEED)?;
    assert!(
        ht_small < HT_NORMAL,
        "Smaller daughter Ht ({ht_small:.4}) should be < feed Ht ({HT_NORMAL:.4})"
    );
    Ok(())
}

/// Daughter hematocrit should always remain in [0, 1] regardless of inputs.
#[test]
fn test_plasma_skimming_bounded() -> cfd_core::error::Result<()> {
    let test_cases = [
        (0.45, 0.1, 30.0, 100.0),
        (0.45, 0.9, 90.0, 100.0),
        (0.80, 0.3, 50.0, 100.0),
        (0.10, 0.5, 100.0, 100.0),
        (0.45, 0.01, 20.0, 200.0),
        (0.45, 0.99, 200.0, 200.0),
    ];

    for &(ht, fq, dd, df) in &test_cases {
        let result = plasma_skimming_hematocrit(ht, fq, diameter_um(dd), diameter_um(df))?;
        assert!(
            (0.0..=1.0).contains(&result),
            "Ht={ht}, fq={fq}, dd={dd}, df={df} → result {result} out of [0,1]"
        );
    }
    Ok(())
}

/// Zero flow fraction means no blood enters this daughter branch,
/// so its hematocrit should be zero.
#[test]
fn test_plasma_skimming_zero_flow_zero_hematocrit() -> cfd_core::error::Result<()> {
    let ht = plasma_skimming_hematocrit(HT_NORMAL, 0.0, diameter_um(50.0), D_FEED)?;
    assert!(
        ht.abs() < 1e-15,
        "Zero flow fraction should give zero hematocrit, got {ht:.10}"
    );
    Ok(())
}

/// Higher feed hematocrit should produce higher daughter hematocrit
/// at the same flow fraction and geometry.
#[test]
fn test_plasma_skimming_increases_with_feed_hematocrit() -> cfd_core::error::Result<()> {
    let fq = 0.4;
    let dd = 70.0;

    let ht_low = plasma_skimming_hematocrit(0.20, fq, diameter_um(dd), D_FEED)?;
    let ht_mid = plasma_skimming_hematocrit(0.35, fq, diameter_um(dd), D_FEED)?;
    let ht_high = plasma_skimming_hematocrit(0.50, fq, diameter_um(dd), D_FEED)?;

    assert!(
        ht_low < ht_mid && ht_mid < ht_high,
        "Daughter Ht should increase with feed Ht: {ht_low:.4} < {ht_mid:.4} < {ht_high:.4}"
    );
    Ok(())
}

#[test]
fn test_pries_phase_separation_enforces_x0_threshold() -> cfd_core::error::Result<()> {
    let result = pries_phase_separation(
        HT_NORMAL,
        0.01,
        diameter_um(30.0),
        diameter_um(90.0),
        diameter_um(20.0),
    )?;
    assert!(
        result.x0 > 0.01,
        "Expected meaningful X0, got {:.4}",
        result.x0
    );
    assert_eq!(result.cell_fraction, 0.0);
    assert_eq!(result.daughter_hematocrit, 0.0);
    Ok(())
}

#[test]
fn test_pries_phase_separation_biases_wider_daughter() -> cfd_core::error::Result<()> {
    let wide = pries_phase_separation(
        HT_NORMAL,
        0.55,
        diameter_um(90.0),
        diameter_um(40.0),
        diameter_um(100.0),
    )?;
    let narrow = pries_phase_separation(
        HT_NORMAL,
        0.45,
        diameter_um(40.0),
        diameter_um(90.0),
        diameter_um(100.0),
    )?;
    assert!(
        wide.cell_fraction > narrow.cell_fraction,
        "Wider daughter should receive more RBC flux than the narrower sibling: wide={:.4}, narrow={:.4}",
        wide.cell_fraction,
        narrow.cell_fraction
    );
    Ok(())
}

#[test]
fn test_checked_pries_matches_legacy_nominal_case() -> cfd_core::error::Result<()> {
    let legacy = pries_phase_separation(
        HT_NORMAL,
        0.55,
        diameter_um(90.0),
        diameter_um(40.0),
        diameter_um(100.0),
    )?;
    let checked = checked_pries_phase_separation(
        HT_NORMAL,
        0.55,
        diameter_um(90.0),
        diameter_um(40.0),
        diameter_um(100.0),
    )?;

    assert!((legacy.cell_fraction - checked.cell_fraction).abs() < 1e-12);
    assert!((legacy.daughter_hematocrit - checked.daughter_hematocrit).abs() < 1e-12);
    Ok(())
}

#[test]
fn test_checked_compact_wrapper_rejects_nonphysical_diameter() {
    let err = checked_plasma_skimming_hematocrit(HT_NORMAL, 0.5, diameter_um(0.0), D_FEED)
        .expect_err("checked plasma-skimming wrapper must reject zero daughter diameter");
    assert!(err.to_string().contains("diameters"));
}

#[test]
fn test_checked_compact_wrapper_matches_legacy_nominal_case() -> cfd_core::error::Result<()> {
    let legacy = plasma_skimming_hematocrit(HT_NORMAL, 0.4, diameter_um(60.0), D_FEED)?;
    let checked = checked_plasma_skimming_hematocrit(HT_NORMAL, 0.4, diameter_um(60.0), D_FEED)?;

    assert!((legacy - checked).abs() < 1e-12);
    Ok(())
}

#[test]
fn test_compact_wrapper_preserves_pries_x0_threshold() -> cfd_core::error::Result<()> {
    let ht = plasma_skimming_hematocrit(HT_NORMAL, 0.01, diameter_um(30.0), diameter_um(20.0))?;
    assert_eq!(
        ht, 0.0,
        "Compact wrapper must inherit the Pries cell-entry threshold"
    );
    Ok(())
}

#[test]
fn test_legacy_wrapper_monotone_for_validation_triplet() -> cfd_core::error::Result<()> {
    let small = plasma_skimming_hematocrit(HT_NORMAL, 0.2, diameter_um(30.0), D_FEED)?;
    let medium = plasma_skimming_hematocrit(HT_NORMAL, 0.4, diameter_um(60.0), D_FEED)?;
    let large = plasma_skimming_hematocrit(HT_NORMAL, 0.6, diameter_um(90.0), D_FEED)?;
    assert!(
        small < medium && medium < large,
        "Legacy wrapper should preserve monotone hematocrit ranking for typical design-space scans: small={small:.4}, medium={medium:.4}, large={large:.4}"
    );
    Ok(())
}
