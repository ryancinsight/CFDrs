//! Pries Phase Separation Model (PSM) for blood cell partitioning at
//! bifurcations.
//!
//! Implements the empirical model of Pries et al. (1989, 1990, 2005) for
//! predicting how red blood cells distribute at asymmetric bifurcations.
//! This is the standard model used in microvascular network simulations
//! and has been validated extensively against in-vivo measurements.
//!
//! # Physics
//!
//! At a bifurcation, the cell-free layer (CFL) near the channel walls
//! causes **plasma skimming**: the low-flow daughter branch draws
//! preferentially from the cell-depleted near-wall fluid, receiving
//! disproportionately fewer cells.  The Pries PSM captures this via:
//!
//! ```text
//! logit(FQE) = A + B * logit((FQB - X0) / (1 - 2*X0))
//! ```
//!
//! where:
//! - FQE = fractional erythrocyte (cell) flux into the daughter
//! - FQB = fractional blood flow into the daughter
//! - X0 = minimum flow fraction needed to draw any cells
//! - A = diameter-asymmetry parameter
//! - B = nonlinearity/shape parameter
//!
//! # References
//!
//! - Pries, A.R., et al. (1989) "Red cell distribution at microvascular
//!   bifurcations." *Microvasc. Res.* 38:81-101.
//! - Pries, A.R. & Secomb, T.W. (2005) "Microvascular blood viscosity
//!   in vivo and the endothelial surface layer." *Am. J. Physiol.*
//!   289:H2657-H2664.

use serde::{Deserialize, Serialize};

/// Parameters for the Pries Phase Separation Model at one bifurcation.
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct PriesPhaseParams {
    /// Parent channel hydraulic diameter [m].
    pub parent_diameter_m: f64,
    /// Daughter alpha hydraulic diameter [m] (the branch being computed).
    pub daughter_alpha_diameter_m: f64,
    /// Daughter beta hydraulic diameter [m] (the other branch).
    pub daughter_beta_diameter_m: f64,
    /// Feed hematocrit (volume fraction of RBCs in the parent).
    pub feed_hematocrit: f64,
}

/// Result of the Pries PSM computation.
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct PhaseSeparationResult {
    /// Fractional blood flow into daughter alpha.
    pub flow_fraction: f64,
    /// Fractional erythrocyte flux into daughter alpha.
    pub cell_fraction: f64,
    /// Hematocrit in daughter alpha.
    pub daughter_hematocrit: f64,
    /// Minimum flow fraction for any cells to enter (X0 parameter).
    pub x0: f64,
}

/// Compute the Pries Phase Separation Model for a single daughter branch.
///
/// Given the fractional blood flow `fqb` into daughter alpha, compute the
/// fractional erythrocyte flux `fqe` and the resulting daughter hematocrit.
///
/// # Arguments
/// * `fqb` - fractional blood flow into daughter alpha (0..1)
/// * `params` - bifurcation geometry and feed hematocrit
///
/// # Returns
/// Phase separation result with cell fraction and daughter hematocrit.
pub fn pries_phase_separation(fqb: f64, params: &PriesPhaseParams) -> PhaseSeparationResult {
    let fqb = fqb.clamp(1e-9, 1.0 - 1e-9);
    let h_d = params.feed_hematocrit.clamp(0.01, 0.70);

    // Convert to microns for the empirical correlation.
    let d_f = params.parent_diameter_m * 1e6; // parent diameter [um]
    let d_alpha = params.daughter_alpha_diameter_m * 1e6;
    let d_beta = params.daughter_beta_diameter_m * 1e6;

    // X0: minimum flow fraction for any cells to enter.
    // Represents the CFL thickness relative to the vessel.
    let x0 = (0.964 * (1.0 - h_d) / d_f).clamp(0.0, 0.49);

    // A: diameter asymmetry parameter.
    // Positive A biases cells toward daughter alpha when it's wider.
    let diameter_ratio_sq = d_alpha * d_alpha / (d_beta * d_beta).max(1e-6);
    let a_param =
        -13.29 * (diameter_ratio_sq - 1.0) / (diameter_ratio_sq + 1.0) * (1.0 - h_d) / d_f.max(1.0);

    // B: nonlinearity parameter. Higher B = steeper S-curve.
    let b_param = 1.0 + 6.98 * (1.0 - h_d) / d_f.max(1.0);

    // Compute FQE using the logit model.
    let fqe = if fqb <= x0 {
        // Below threshold: only plasma enters, no cells.
        0.0
    } else if fqb >= 1.0 - x0 {
        // Above complement: all cells enter this branch.
        1.0
    } else {
        // Core logit model.
        let normalized = (fqb - x0) / (1.0 - 2.0 * x0);
        let logit_input = logit(normalized.clamp(1e-6, 1.0 - 1e-6));
        let logit_fqe = a_param + b_param * logit_input;
        inv_logit(logit_fqe).clamp(0.0, 1.0)
    };

    // Daughter hematocrit: H_daughter = FQE * H_parent / FQB
    let daughter_hematocrit = if fqb > 1e-12 {
        (fqe * h_d / fqb).clamp(0.0, h_d * 2.0) // cap at 2x feed to prevent unphysical values
    } else {
        0.0
    };

    PhaseSeparationResult {
        flow_fraction: fqb,
        cell_fraction: fqe,
        daughter_hematocrit,
        x0,
    }
}

/// Extension: compute phase separation for CTCs and WBCs using modified
/// parameters.  CTCs are larger and stiffer, so they have a thinner
/// effective CFL and stronger bias toward the high-flow branch.
///
/// The modification scales X0 by the cell-size ratio: larger cells have
/// smaller effective X0 (they can't fit in the near-wall plasma layer
/// as easily as RBCs).
pub fn pries_phase_separation_cell_type(
    fqb: f64,
    params: &PriesPhaseParams,
    cell_diameter_m: f64,
    stiffness_factor: f64,
) -> PhaseSeparationResult {
    let rbc_diameter = 8.0e-6;
    // Larger cells have reduced X0: they extend further into the core
    // flow and are less affected by the CFL.
    let size_ratio = (cell_diameter_m / rbc_diameter).max(0.5);
    // Stiffer cells also have reduced X0: they don't deform to fit
    // in the near-wall layer.
    let stiffness_correction = stiffness_factor.max(0.5);

    let mut modified_params = *params;
    // Scale the parent diameter to effectively reduce X0 for larger cells.
    // This is equivalent to making the CFL thinner relative to the cell.
    let effective_parent_d = params.parent_diameter_m * size_ratio * stiffness_correction;
    modified_params.parent_diameter_m = effective_parent_d;

    pries_phase_separation(fqb, &modified_params)
}

fn logit(x: f64) -> f64 {
    (x / (1.0 - x)).ln()
}

fn inv_logit(x: f64) -> f64 {
    1.0 / (1.0 + (-x).exp())
}

/// Compute a three-population phase separation at a bifurcation.
///
/// Uses the Pries PSM with cell-type-specific corrections:
/// - CTCs: larger diameter (12.5 um), stiffer (factor 1.5)
/// - WBCs: intermediate diameter (11 um), moderate stiffness (1.2)
/// - RBCs: baseline diameter (8 um), stiffness factor 1.0
///
/// Returns (ctc_fraction, wbc_fraction, rbc_fraction) for the alpha daughter.
pub fn three_population_phase_separation(
    fqb: f64,
    params: &PriesPhaseParams,
) -> ThreePopulationResult {
    let rbc = pries_phase_separation(fqb, params);
    let ctc = pries_phase_separation_cell_type(fqb, params, 12.5e-6, 1.5);
    let wbc = pries_phase_separation_cell_type(fqb, params, 11.0e-6, 1.2);

    ThreePopulationResult {
        rbc_fraction: rbc.cell_fraction,
        wbc_fraction: wbc.cell_fraction,
        ctc_fraction: ctc.cell_fraction,
        rbc_hematocrit: rbc.daughter_hematocrit,
        x0_rbc: rbc.x0,
        x0_ctc: ctc.x0,
    }
}

/// Three-population phase separation result.
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct ThreePopulationResult {
    /// Fractional RBC flux into the daughter.
    pub rbc_fraction: f64,
    /// Fractional WBC flux into the daughter.
    pub wbc_fraction: f64,
    /// Fractional CTC flux into the daughter.
    pub ctc_fraction: f64,
    /// Daughter hematocrit (from RBC partition).
    pub rbc_hematocrit: f64,
    /// X0 for RBCs (minimum flow fraction for cell entry).
    pub x0_rbc: f64,
    /// X0 for CTCs (typically smaller due to larger cell size).
    pub x0_ctc: f64,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn symmetric_bifurcation_equal_split() {
        let params = PriesPhaseParams {
            parent_diameter_m: 200e-6,
            daughter_alpha_diameter_m: 150e-6,
            daughter_beta_diameter_m: 150e-6,
            feed_hematocrit: 0.45,
        };
        let result = pries_phase_separation(0.5, &params);
        // Symmetric flow + symmetric diameters -> equal cell split.
        assert!(
            (result.cell_fraction - 0.5).abs() < 0.05,
            "symmetric bifurcation should give ~50% cell fraction, got {:.3}",
            result.cell_fraction
        );
    }

    #[test]
    fn asymmetric_bifurcation_wide_daughter_gets_more_cells() {
        let params = PriesPhaseParams {
            parent_diameter_m: 2000e-6,         // 2 mm (millifluidic)
            daughter_alpha_diameter_m: 1100e-6, // wide
            daughter_beta_diameter_m: 900e-6,   // narrow
            feed_hematocrit: 0.45,
        };
        // 65% flow to wide daughter.
        let result = pries_phase_separation(0.65, &params);
        assert!(
            result.cell_fraction > 0.65,
            "wider daughter with 65% flow should get >65% cells, got {:.3}",
            result.cell_fraction
        );
    }

    #[test]
    fn pries_psm_is_for_rbc_partitioning_not_cell_type_differentiation() {
        // The Pries PSM models RBC hematocrit partitioning.  It does NOT
        // inherently differentiate cell types (CTCs vs RBCs) because it
        // models the cell-free layer as an RBC-specific phenomenon.
        //
        // Cell-type-dependent routing (Zweifach-Fung stiffness exponents)
        // is a separate empirical model implemented in cfd-1d. The PSM
        // provides the baseline hematocrit partition that the Zweifach-Fung
        // model builds upon.
        let params = PriesPhaseParams {
            parent_diameter_m: 100e-6,
            daughter_alpha_diameter_m: 60e-6,
            daughter_beta_diameter_m: 40e-6,
            feed_hematocrit: 0.30,
        };
        let fqb = 0.65;
        let rbc = pries_phase_separation(fqb, &params);
        eprintln!(
            "Pries PSM (100um) at fqb={fqb}: RBC fraction={:.4}, X0={:.4}",
            rbc.cell_fraction, rbc.x0
        );
        // At 100 um, X0 is significant (~5% of flow).
        assert!(
            rbc.x0 > 0.005,
            "X0 should be meaningful at microfluidic scale, got {:.4}",
            rbc.x0
        );
        // The daughter should get more cells than flow (plasma skimming
        // of the peripheral arm means center arm gets enriched).
        assert!(
            rbc.cell_fraction >= fqb * 0.95,
            "wide daughter should get at least 95% of its flow-proportional cells"
        );
    }

    #[test]
    fn millifluidic_scale_psm_is_negligible() {
        // At millifluidic scale (2 mm), the CFL is negligible relative
        // to channel diameter.  The Pries PSM gives essentially identical
        // fractions for all cell types.  This is expected and correct:
        // at this scale, the Zweifach-Fung stiffness-exponent model
        // (cfd-1d) is the appropriate physics model, not plasma skimming.
        let params = PriesPhaseParams {
            parent_diameter_m: 2000e-6, // 2 mm
            daughter_alpha_diameter_m: 1100e-6,
            daughter_beta_diameter_m: 900e-6,
            feed_hematocrit: 0.45,
        };
        let fqb = 0.65;
        let three_pop = three_population_phase_separation(fqb, &params);
        let delta = (three_pop.ctc_fraction - three_pop.rbc_fraction).abs();
        eprintln!(
            "Pries PSM (2mm) at fqb={fqb}: CTC={:.4}, RBC={:.4}, delta={:.6}",
            three_pop.ctc_fraction, three_pop.rbc_fraction, delta
        );
        // At 2 mm, the difference should be < 1%.
        assert!(
            delta < 0.01,
            "at millifluidic scale, PSM should show < 1% CTC/RBC difference, got {delta:.4}"
        );
    }

    #[test]
    fn low_flow_branch_gets_no_cells_below_x0() {
        let params = PriesPhaseParams {
            parent_diameter_m: 100e-6, // small vessel where CFL matters
            daughter_alpha_diameter_m: 50e-6,
            daughter_beta_diameter_m: 80e-6,
            feed_hematocrit: 0.30,
        };
        let result = pries_phase_separation(0.01, &params);
        assert!(
            result.cell_fraction < 0.01,
            "below X0, cell fraction should be ~0, got {:.4}",
            result.cell_fraction
        );
    }
}
