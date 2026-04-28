//! Cross-track Milestone 12 therapy utility.
//!
//! Track-local objective scores are valid only within their search modes:
//! Option 1 optimizes acoustic residence/separation, Option 2 optimizes
//! venturi cavitation selectivity, and GA optimizes architecture-preserving
//! Dean/venturi refinement. This module provides the report-level score used
//! when those tracks are compared in one figure or table.
//!
//! # Theorem
//! The cross-track utility is comparable across Option 1, Option 2, and GA
//! outputs because every design is projected onto the same bounded physical
//! coordinates: cancer-selective cavitation delivery, three-population routing,
//! healthy-cell protection, treatment dwell, pediatric ECV margin, and hard
//! safety compliance.
//!
//! **Proof sketch**
//! Each component is a dimensionless value in `[0, 1]` derived from the same
//! [`crate::SdtMetrics`] fields regardless of optimization track. Ultrasound
//! SDT receives acoustic-cavitation delivery credit only when the treatment mode
//! contains ultrasound actuation and the channel has non-zero half-wavelength
//! resonance. Hydrodynamic SDT receives delivery credit only when finite
//! cavitating venturi geometry is active. Designs with neither mechanism are
//! mapped through the explicitly capped separation-support branch. Cavitating
//! designs use one weighted sum plus one geometric-mean coupling term over the
//! same component vector, then clamp to `[0, 1]`. Therefore two designs with the
//! same component values receive the same utility independent of their
//! generating objective, while ultrasound-only designs are not misclassified as
//! non-cavitating merely because they have no venturi throat.

use crate::constraints::{PEDIATRIC_BLOOD_VOLUME_ML_PER_KG, PEDIATRIC_REFERENCE_WEIGHT_KG};
use crate::SdtMetrics;

const RESIDENCE_REFERENCE_S: f64 = 0.020;
const SEPARATION_SUPPORT_CAP: f64 = 0.35;

/// Component-level explanation of the cross-track therapy utility.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TherapyUtilityComponents {
    pub cavitation_delivery: f64,
    pub routing: f64,
    pub healthy_cell_protection: f64,
    pub residence: f64,
    pub pediatric_ecv_margin: f64,
    pub safety: f64,
    pub synergy: f64,
    pub separation_support_only: bool,
}

/// Compute the report-level Milestone 12 therapy utility.
#[must_use]
pub fn milestone12_therapy_utility(metrics: &SdtMetrics) -> f64 {
    let components = therapy_utility_components(metrics);
    if components.separation_support_only {
        return (SEPARATION_SUPPORT_CAP
            * separation_support(
                components.routing,
                components.healthy_cell_protection,
                components.residence,
                components.safety,
            ))
        .clamp(0.0, SEPARATION_SUPPORT_CAP);
    }

    let base = 0.30 * components.cavitation_delivery
        + 0.25 * components.healthy_cell_protection
        + 0.20 * components.routing
        + 0.15 * components.residence
        + 0.05 * components.pediatric_ecv_margin
        + 0.05 * components.safety;

    (base + components.synergy).clamp(0.0, 1.0)
}

/// Return the normalized physical components used by
/// [`milestone12_therapy_utility`].
#[must_use]
pub fn therapy_utility_components(metrics: &SdtMetrics) -> TherapyUtilityComponents {
    let hydrodynamic_delivery = hydrodynamic_cavitation_delivery(metrics);
    let acoustic_delivery = acoustic_cavitation_delivery(metrics);
    let cavitation_delivery = hydrodynamic_delivery.max(acoustic_delivery);
    let has_cavitation_delivery = cavitation_delivery > 0.0;

    let routing = metrics.three_pop_sep_efficiency.clamp(0.0, 1.0);
    let healthy_cell_protection = metrics.healthy_cell_protection_index.clamp(0.0, 1.0);
    let residence = (metrics.treatment_zone_dwell_time_s / RESIDENCE_REFERENCE_S).clamp(0.0, 1.0);
    let pediatric_ecv_margin = pediatric_ecv_margin(metrics.total_ecv_ml);
    let safety = if metrics.pressure_feasible
        && metrics.plate_fits
        && metrics.fda_main_compliant
        && metrics.fda_overall_compliant
    {
        1.0
    } else {
        0.0
    };
    let synergy = if has_cavitation_delivery {
        0.10 * (cavitation_delivery * healthy_cell_protection * routing * residence).powf(0.25)
    } else {
        0.0
    };

    TherapyUtilityComponents {
        cavitation_delivery,
        routing,
        healthy_cell_protection,
        residence,
        pediatric_ecv_margin,
        safety,
        synergy,
        separation_support_only: !has_cavitation_delivery,
    }
}

fn hydrodynamic_cavitation_delivery(metrics: &SdtMetrics) -> f64 {
    let has_cavitating_venturi = metrics.active_venturi_throat_count > 0
        && metrics.cavitation_number.is_finite()
        && metrics.cavitation_number < 1.0
        && metrics.cancer_targeted_cavitation > 0.0;
    if has_cavitating_venturi {
        (0.45 * metrics.cancer_targeted_cavitation.clamp(0.0, 1.0)
            + 0.35 * metrics.serial_cavitation_dose_fraction.clamp(0.0, 1.0)
            + 0.20 * metrics.cavitation_intensity.clamp(0.0, 1.0))
        .clamp(0.0, 1.0)
    } else {
        0.0
    }
}

fn acoustic_cavitation_delivery(metrics: &SdtMetrics) -> f64 {
    if !metrics.treatment_zone_mode.contains("Ultrasound") {
        return 0.0;
    }
    (metrics.acoustic_resonance_factor.clamp(0.0, 1.0)
        * metrics.cancer_center_fraction.clamp(0.0, 1.0)
        * metrics.serial_cavitation_dose_fraction.clamp(0.0, 1.0))
    .clamp(0.0, 1.0)
}

fn separation_support(
    routing: f64,
    healthy_cell_protection: f64,
    residence: f64,
    safety: f64,
) -> f64 {
    (0.35 * routing + 0.35 * healthy_cell_protection + 0.20 * residence + 0.10 * safety)
        .clamp(0.0, 1.0)
}

fn pediatric_ecv_margin(ecv_ml: f64) -> f64 {
    let limit_ml = PEDIATRIC_REFERENCE_WEIGHT_KG * PEDIATRIC_BLOOD_VOLUME_ML_PER_KG * 0.10;
    (1.0 - ecv_ml.max(0.0) / limit_ml.max(1.0e-12)).clamp(0.0, 1.0)
}

#[cfg(test)]
mod tests {
    use super::{milestone12_therapy_utility, therapy_utility_components, SEPARATION_SUPPORT_CAP};
    use crate::domain::fixtures::{canonical_option2_candidate, operating_point};
    use crate::reporting::compute_blueprint_report_metrics;

    #[test]
    fn non_cavitating_design_is_capped_as_separation_support() {
        let candidate = canonical_option2_candidate(
            "therapy-utility-cap",
            operating_point(1.0e-6, 5_000.0, 0.2),
        );
        let mut metrics = compute_blueprint_report_metrics(&candidate).expect("metrics");
        metrics.active_venturi_throat_count = 0;
        metrics.cavitation_number = f64::NAN;
        metrics.treatment_zone_mode = "SeparationOnly".to_string();
        metrics.acoustic_resonance_factor = 0.0;
        metrics.cancer_targeted_cavitation = 0.0;
        metrics.serial_cavitation_dose_fraction = 0.0;
        metrics.cavitation_intensity = 0.0;

        let components = therapy_utility_components(&metrics);
        let utility = milestone12_therapy_utility(&metrics);

        assert!(components.separation_support_only);
        assert_eq!(components.cavitation_delivery, 0.0);
        assert!(
            utility <= SEPARATION_SUPPORT_CAP,
            "non-cavitating separation-only utility {utility} must not exceed cap {SEPARATION_SUPPORT_CAP}"
        );
    }

    #[test]
    fn ultrasound_only_design_receives_acoustic_cavitation_delivery() {
        let candidate = crate::domain::fixtures::canonical_option1_candidate(
            "therapy-utility-acoustic",
            operating_point(2.0e-6, 30_000.0, 0.18),
        );
        let mut metrics = compute_blueprint_report_metrics(&candidate).expect("metrics");
        metrics.active_venturi_throat_count = 0;
        metrics.cavitation_number = f64::INFINITY;
        metrics.cavitation_intensity = 0.0;
        metrics.cancer_targeted_cavitation = 0.0;
        metrics.treatment_zone_mode = "UltrasoundOnly".to_string();
        metrics.acoustic_resonance_factor = 0.856_103_896_103_896;
        metrics.cancer_center_fraction = 0.886_010_994_883_551_5;
        metrics.serial_cavitation_dose_fraction = 0.136_125;

        let components = therapy_utility_components(&metrics);
        let utility = milestone12_therapy_utility(&metrics);
        let expected_delivery = metrics.acoustic_resonance_factor
            * metrics.cancer_center_fraction
            * metrics.serial_cavitation_dose_fraction;

        assert!(!components.separation_support_only);
        assert!(
            (components.cavitation_delivery - expected_delivery).abs() <= 1.0e-12,
            "acoustic delivery mismatch: got {}, expected {}",
            components.cavitation_delivery,
            expected_delivery
        );
        assert!(
            utility > SEPARATION_SUPPORT_CAP,
            "ultrasound-cavitating utility {utility} must not be capped as separation-only"
        );
    }

    #[test]
    fn utility_increases_with_cavitation_dose_at_fixed_protection() {
        let candidate = canonical_option2_candidate(
            "therapy-utility-cavitation-monotone",
            operating_point(2.0e-6, 300_000.0, 0.18),
        );
        let mut low = compute_blueprint_report_metrics(&candidate).expect("metrics");
        low.active_venturi_throat_count = 1;
        low.cavitation_number = 0.5;
        low.cancer_targeted_cavitation = 0.20;
        low.serial_cavitation_dose_fraction = 0.25;
        low.cavitation_intensity = 0.25;

        let mut high = low.clone();
        high.cancer_targeted_cavitation = 0.60;
        high.serial_cavitation_dose_fraction = 0.75;
        high.cavitation_intensity = 0.75;

        let low_score = milestone12_therapy_utility(&low);
        let high_score = milestone12_therapy_utility(&high);

        assert!(
            high_score > low_score,
            "higher cavitation dose should raise utility: low={low_score}, high={high_score}"
        );
    }
}
