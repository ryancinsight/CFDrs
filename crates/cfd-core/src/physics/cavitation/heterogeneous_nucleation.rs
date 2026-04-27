//! Selective heterogeneous cavitation thresholds for mixed cellular populations.
//!
//! # Theorem — Ordered selective Blake thresholds
//!
//! For a physically admissible population set with positive radius, positive
//! ambient overpressure above vapor pressure, nonnegative interfacial tension,
//! and bounded nuclei weighting, each population admits a finite effective
//! inception threshold. If one population has lower effective interfacial
//! tension, lower membrane stiffness, or larger effective seed radius than
//! another while holding the remaining terms fixed, its effective inception
//! threshold is greater and therefore cavitation occurs earlier under the same
//! external tensile excursion.
//!
//! **Proof sketch**: the classical Blake threshold rises toward vapor pressure
//! as radius increases or surface tension decreases. The model below maps
//! membrane compliance into a bounded softening of effective tension and a
//! bounded increase in effective seed radius, then adds a bounded nuclei-driven
//! pressure lift. Each transformation is monotone in the physically intended
//! direction, so the ordering of effective thresholds is preserved.

use crate::error::{Error, Result};
use crate::physics::cavitation::nuclei_transport::NUCLEI_VAPOR_PRESSURE_BOOST_PA_PER_UNIT_FRACTION;
use nalgebra::RealField;
use serde::{Deserialize, Serialize};

const REFERENCE_MEMBRANE_STIFFNESS_PA: f64 = 100_000.0;
const MAX_COMPLIANCE_GAIN: f64 = 4.0;
const MIN_COMPLIANCE_GAIN: f64 = 0.25;
const MAX_RADIUS_AMPLIFICATION: f64 = 0.35;
const MAX_SEED_DENSITY_FACTOR: f64 = 0.15;

/// Canonical identity for selective cavitation screening.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum CellPopulationIdentity {
    /// Healthy red blood cells.
    HealthyRbc,
    /// Healthy white blood cells.
    HealthyWbc,
    /// Circulating tumor cells or analogous target cells.
    CirculatingTumorCell,
    /// A generic targeted population.
    GenericTarget,
    /// A generic healthy or non-target population.
    GenericHealthy,
}

impl CellPopulationIdentity {
    /// Whether this population is a therapeutic target.
    #[must_use]
    pub const fn is_target(self) -> bool {
        matches!(self, Self::CirculatingTumorCell | Self::GenericTarget)
    }
}

/// Mechanical state entering the selective cavitation model.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct CellMechanicalState {
    /// Effective membrane stiffness [Pa].
    pub membrane_stiffness_pa: f64,
    /// Effective interfacial tension [N/m].
    pub interfacial_tension_n_m: f64,
    /// Characteristic nuclei seed radius [m].
    pub particle_radius_m: f64,
    /// Optional bounded deformability multiplier; 1.0 means neutral.
    #[serde(default = "unit_f64")]
    pub deformability_factor: f64,
}

/// Nucleation state entering the selective cavitation model.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct PopulationNucleationState {
    /// Mixture fraction of the population.
    pub volume_fraction: f64,
    /// Upstream nuclei fraction [0, 1].
    #[serde(default)]
    pub upstream_nuclei_fraction: f64,
    /// Seed density proxy [0, +inf), internally bounded.
    #[serde(default)]
    pub seed_density_factor: f64,
    /// Relative inception weighting [0, 1].
    #[serde(default = "unit_f64")]
    pub inception_weight: f64,
}

/// Population definition for selective cavitation screening.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct SelectiveCavitationPopulation {
    /// Population identity.
    pub identity: CellPopulationIdentity,
    /// Human-readable label for reporting.
    pub label: String,
    /// Mechanical state.
    pub mechanical_state: CellMechanicalState,
    /// Nucleation state.
    pub nucleation_state: PopulationNucleationState,
}

/// Input state for selective cavitation evaluation.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct SelectiveCavitationInput {
    /// Base vapor pressure [Pa].
    pub base_vapor_pressure_pa: f64,
    /// Ambient static pressure [Pa].
    pub ambient_pressure_pa: f64,
    /// Fluid density [kg/m^3].
    pub density_kg_m3: f64,
    /// Population set.
    pub populations: Vec<SelectiveCavitationPopulation>,
}

/// Per-population threshold report.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct PopulationCavitationThreshold {
    /// Population identity.
    pub identity: CellPopulationIdentity,
    /// Human-readable label.
    pub label: String,
    /// Compliance gain used internally.
    pub compliance_gain: f64,
    /// Effective seed radius [m].
    pub effective_particle_radius_m: f64,
    /// Effective interfacial tension [N/m].
    pub effective_interfacial_tension_n_m: f64,
    /// Classical Blake threshold before nuclei lifting [Pa].
    pub blake_threshold_pressure_pa: f64,
    /// Effective threshold after nuclei lifting [Pa].
    pub effective_threshold_pressure_pa: f64,
    /// Weight used in mixture aggregation.
    pub mixture_weight: f64,
}

/// Aggregate selective cavitation result.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct SelectiveCavitationResult {
    /// Per-population thresholds sorted from highest to lowest effective threshold.
    pub population_thresholds: Vec<PopulationCavitationThreshold>,
    /// Volume-weighted mixture inception threshold [Pa].
    pub mixture_inception_threshold_pa: f64,
    /// Highest effective threshold among target populations [Pa].
    pub target_inception_threshold_pa: f64,
    /// Highest effective threshold among healthy populations [Pa].
    pub healthy_inception_threshold_pa: f64,
    /// Positive means the best target cavitates earlier than the best healthy population.
    pub selectivity_margin_pa: f64,
    /// Dominant selective population if one exists.
    pub dominant_selective_population: Option<CellPopulationIdentity>,
    /// Dominant selective label if one exists.
    pub dominant_selective_label: Option<String>,
}

/// Legacy coarse population type retained as a compatibility adapter.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct CellularPopulation<T: RealField + Copy> {
    /// Volume fraction of instances.
    pub volume_fraction: T,
    /// Absolute membrane stiffness scaling [Pa].
    pub membrane_stiffness_pa: T,
    /// Interfacial tension of the cellular envelope [N/m].
    pub interfacial_tension_n_m: T,
    /// Nominal particle / cell radius acting as initial nuclei seed [m].
    pub particle_radius_m: T,
}

fn unit_f64() -> f64 {
    1.0
}

fn compliance_gain(stiffness_pa: f64, deformability_factor: f64) -> f64 {
    let stiffness_term = (REFERENCE_MEMBRANE_STIFFNESS_PA / stiffness_pa.max(1.0))
        .sqrt()
        .clamp(MIN_COMPLIANCE_GAIN, MAX_COMPLIANCE_GAIN);
    (stiffness_term * deformability_factor.max(0.25))
        .clamp(MIN_COMPLIANCE_GAIN, MAX_COMPLIANCE_GAIN)
}

fn blake_threshold_pressure_pa(
    vapor_pressure_pa: f64,
    ambient_pressure_pa: f64,
    sigma_n_m: f64,
    seed_radius_m: f64,
) -> f64 {
    if sigma_n_m <= 0.0 {
        return vapor_pressure_pa;
    }

    let gas_overpressure =
        (ambient_pressure_pa - vapor_pressure_pa + 2.0 * sigma_n_m / seed_radius_m).max(0.0);
    let critical_radius = ((3.0 * seed_radius_m.powi(3) * gas_overpressure) / (2.0 * sigma_n_m))
        .sqrt()
        .max(seed_radius_m);
    // Classical Blake threshold:
    // P_B = P_v - 4σ / (3 R_B)
    vapor_pressure_pa - (4.0 * sigma_n_m) / (3.0 * critical_radius)
}

/// Validate a selective cavitation input contract.
pub fn validate_selective_cavitation_input(input: &SelectiveCavitationInput) -> Result<()> {
    if !input.base_vapor_pressure_pa.is_finite()
        || !input.ambient_pressure_pa.is_finite()
        || !input.density_kg_m3.is_finite()
    {
        return Err(Error::InvalidConfiguration(
            "selective cavitation scalar inputs must be finite".to_string(),
        ));
    }
    if input.base_vapor_pressure_pa < 0.0 {
        return Err(Error::InvalidConfiguration(
            "base vapor pressure must be nonnegative".to_string(),
        ));
    }
    if input.ambient_pressure_pa <= input.base_vapor_pressure_pa {
        return Err(Error::InvalidConfiguration(
            "ambient pressure must exceed vapor pressure for Blake-threshold evaluation"
                .to_string(),
        ));
    }
    if input.density_kg_m3 <= 0.0 {
        return Err(Error::InvalidConfiguration(
            "fluid density must be positive".to_string(),
        ));
    }
    if input.populations.is_empty() {
        return Err(Error::InvalidConfiguration(
            "selective cavitation requires at least one population".to_string(),
        ));
    }

    for population in &input.populations {
        let mechanics = population.mechanical_state;
        let nucleation = population.nucleation_state;
        if population.label.trim().is_empty() {
            return Err(Error::InvalidConfiguration(
                "population labels must be non-empty".to_string(),
            ));
        }
        if !mechanics.membrane_stiffness_pa.is_finite()
            || !mechanics.interfacial_tension_n_m.is_finite()
            || !mechanics.particle_radius_m.is_finite()
            || !mechanics.deformability_factor.is_finite()
            || !nucleation.volume_fraction.is_finite()
            || !nucleation.upstream_nuclei_fraction.is_finite()
            || !nucleation.seed_density_factor.is_finite()
            || !nucleation.inception_weight.is_finite()
        {
            return Err(Error::InvalidConfiguration(
                "population cavitation inputs must be finite".to_string(),
            ));
        }
        if mechanics.membrane_stiffness_pa <= 0.0
            || mechanics.interfacial_tension_n_m < 0.0
            || mechanics.particle_radius_m <= 0.0
            || mechanics.deformability_factor <= 0.0
        {
            return Err(Error::InvalidConfiguration(
                "population mechanics must be physically positive".to_string(),
            ));
        }
        if !(0.0..=1.0).contains(&nucleation.volume_fraction)
            || !(0.0..=1.0).contains(&nucleation.upstream_nuclei_fraction)
            || !(0.0..=1.0).contains(&nucleation.inception_weight)
            || nucleation.seed_density_factor < 0.0
        {
            return Err(Error::InvalidConfiguration(
                "population nucleation weights must be bounded and nonnegative".to_string(),
            ));
        }
    }

    Ok(())
}

/// Evaluate selective cavitation thresholds and ordering.
pub fn evaluate_selective_cavitation_thresholds(
    input: &SelectiveCavitationInput,
) -> Result<SelectiveCavitationResult> {
    validate_selective_cavitation_input(input)?;

    let mut thresholds = Vec::with_capacity(input.populations.len());
    let mut weighted_sum = 0.0;
    let mut total_weight = 0.0;
    let mut target_best = f64::NEG_INFINITY;
    let mut healthy_best = f64::NEG_INFINITY;
    let mut dominant_identity = None;
    let mut dominant_label = None;

    for population in &input.populations {
        let mechanics = population.mechanical_state;
        let nucleation = population.nucleation_state;
        let gain = compliance_gain(
            mechanics.membrane_stiffness_pa,
            mechanics.deformability_factor,
        );
        let radius_amplification = 1.0 + MAX_RADIUS_AMPLIFICATION * (gain - 1.0).max(0.0);
        let sigma_relaxation = gain.sqrt().max(1.0);
        let effective_radius = mechanics.particle_radius_m * radius_amplification;
        let effective_sigma = mechanics.interfacial_tension_n_m / sigma_relaxation;
        let base_threshold = blake_threshold_pressure_pa(
            input.base_vapor_pressure_pa,
            input.ambient_pressure_pa,
            effective_sigma,
            effective_radius,
        );
        let bounded_seed_density = nucleation.seed_density_factor.min(1.0);
        let nuclei_lift = nucleation.inception_weight
            * (nucleation.upstream_nuclei_fraction
                + MAX_SEED_DENSITY_FACTOR * bounded_seed_density)
            * NUCLEI_VAPOR_PRESSURE_BOOST_PA_PER_UNIT_FRACTION;
        let effective_threshold = base_threshold + nuclei_lift;

        weighted_sum += effective_threshold * nucleation.volume_fraction;
        total_weight += nucleation.volume_fraction;

        if population.identity.is_target() && effective_threshold > target_best {
            target_best = effective_threshold;
        }
        if !population.identity.is_target() && effective_threshold > healthy_best {
            healthy_best = effective_threshold;
        }

        thresholds.push(PopulationCavitationThreshold {
            identity: population.identity,
            label: population.label.clone(),
            compliance_gain: gain,
            effective_particle_radius_m: effective_radius,
            effective_interfacial_tension_n_m: effective_sigma,
            blake_threshold_pressure_pa: base_threshold,
            effective_threshold_pressure_pa: effective_threshold,
            mixture_weight: nucleation.volume_fraction,
        });
    }

    thresholds.sort_by(|lhs, rhs| {
        rhs.effective_threshold_pressure_pa
            .partial_cmp(&lhs.effective_threshold_pressure_pa)
            .unwrap_or(std::cmp::Ordering::Equal)
    });

    if let Some(top) = thresholds.first() {
        dominant_identity = Some(top.identity);
        dominant_label = Some(top.label.clone());
    }

    let mixture_inception_threshold_pa = if total_weight > 0.0 {
        weighted_sum / total_weight
    } else {
        input.base_vapor_pressure_pa
    };
    let target_inception_threshold_pa = if target_best.is_finite() {
        target_best
    } else {
        mixture_inception_threshold_pa
    };
    let healthy_inception_threshold_pa = if healthy_best.is_finite() {
        healthy_best
    } else {
        mixture_inception_threshold_pa
    };
    let selectivity_margin_pa = target_inception_threshold_pa - healthy_inception_threshold_pa;

    Ok(SelectiveCavitationResult {
        population_thresholds: thresholds,
        mixture_inception_threshold_pa,
        target_inception_threshold_pa,
        healthy_inception_threshold_pa,
        selectivity_margin_pa,
        dominant_selective_population: dominant_identity,
        dominant_selective_label: dominant_label,
    })
}

/// Return thresholds sorted from most to least cavitation-prone.
pub fn rank_population_selectivity(
    input: &SelectiveCavitationInput,
) -> Result<Vec<PopulationCavitationThreshold>> {
    Ok(evaluate_selective_cavitation_thresholds(input)?.population_thresholds)
}

/// Compatibility adapter for the legacy coarse mixture helper.
pub fn heterogeneous_inception_threshold_pa<T: RealField + Copy>(
    base_vapor_pressure_pa: T,
    ambient_pressure_pa: T,
    populations: &[CellularPopulation<T>],
) -> T {
    let mapped = populations
        .iter()
        .enumerate()
        .map(|(index, pop)| SelectiveCavitationPopulation {
            identity: CellPopulationIdentity::GenericHealthy,
            label: format!("population_{index}"),
            mechanical_state: CellMechanicalState {
                membrane_stiffness_pa: pop
                    .membrane_stiffness_pa
                    .to_subset()
                    .unwrap_or(REFERENCE_MEMBRANE_STIFFNESS_PA),
                interfacial_tension_n_m: pop.interfacial_tension_n_m.to_subset().unwrap_or(0.05),
                particle_radius_m: pop.particle_radius_m.to_subset().unwrap_or(5.0e-6),
                deformability_factor: 1.0,
            },
            nucleation_state: PopulationNucleationState {
                volume_fraction: pop.volume_fraction.to_subset().unwrap_or(0.0),
                upstream_nuclei_fraction: 0.0,
                seed_density_factor: 0.0,
                inception_weight: 1.0,
            },
        })
        .collect::<Vec<_>>();

    let input = SelectiveCavitationInput {
        base_vapor_pressure_pa: base_vapor_pressure_pa.to_subset().unwrap_or(0.0),
        ambient_pressure_pa: ambient_pressure_pa.to_subset().unwrap_or(101_325.0),
        density_kg_m3: 1_000.0,
        populations: mapped,
    };

    let threshold = evaluate_selective_cavitation_thresholds(&input).map_or_else(
        |_| input.base_vapor_pressure_pa,
        |result| result.mixture_inception_threshold_pa,
    );
    T::from_f64(threshold).unwrap_or(base_vapor_pressure_pa)
}

#[cfg(test)]
mod tests {
    use super::*;
    use proptest::prelude::*;

    fn healthy_population() -> SelectiveCavitationPopulation {
        SelectiveCavitationPopulation {
            identity: CellPopulationIdentity::HealthyRbc,
            label: "healthy_rbc".to_string(),
            mechanical_state: CellMechanicalState {
                membrane_stiffness_pa: 120_000.0,
                interfacial_tension_n_m: 0.07,
                particle_radius_m: 4.0e-6,
                deformability_factor: 1.0,
            },
            nucleation_state: PopulationNucleationState {
                volume_fraction: 0.8,
                upstream_nuclei_fraction: 0.01,
                seed_density_factor: 0.1,
                inception_weight: 1.0,
            },
        }
    }

    fn target_population() -> SelectiveCavitationPopulation {
        SelectiveCavitationPopulation {
            identity: CellPopulationIdentity::CirculatingTumorCell,
            label: "ctc".to_string(),
            mechanical_state: CellMechanicalState {
                membrane_stiffness_pa: 20_000.0,
                interfacial_tension_n_m: 0.03,
                particle_radius_m: 9.0e-6,
                deformability_factor: 1.25,
            },
            nucleation_state: PopulationNucleationState {
                volume_fraction: 0.2,
                upstream_nuclei_fraction: 0.04,
                seed_density_factor: 0.35,
                inception_weight: 1.0,
            },
        }
    }

    #[test]
    fn selective_target_cavitates_earlier_than_healthy_population() {
        let result = evaluate_selective_cavitation_thresholds(&SelectiveCavitationInput {
            base_vapor_pressure_pa: 3_170.0,
            ambient_pressure_pa: 101_325.0,
            density_kg_m3: 1_000.0,
            populations: vec![healthy_population(), target_population()],
        })
        .expect("selective cavitation should evaluate");

        assert!(result.selectivity_margin_pa > 0.0);
        assert_eq!(
            result.dominant_selective_population,
            Some(CellPopulationIdentity::CirculatingTumorCell)
        );
    }

    #[test]
    fn blake_threshold_matches_classic_four_thirds_factor() {
        let vapor_pressure_pa = 3_170.0_f64;
        let ambient_pressure_pa = 101_325.0_f64;
        let sigma_n_m = 0.072_f64;
        let seed_radius_m = 5.0e-6_f64;

        let gas_overpressure = (ambient_pressure_pa - vapor_pressure_pa
            + 2.0_f64 * sigma_n_m / seed_radius_m)
            .max(0.0_f64);
        let critical_radius = ((3.0_f64 * seed_radius_m.powi(3) * gas_overpressure)
            / (2.0_f64 * sigma_n_m))
            .sqrt()
            .max(seed_radius_m);
        let expected = vapor_pressure_pa - (4.0_f64 * sigma_n_m) / (3.0_f64 * critical_radius);

        let actual = blake_threshold_pressure_pa(
            vapor_pressure_pa,
            ambient_pressure_pa,
            sigma_n_m,
            seed_radius_m,
        );

        assert!((actual - expected).abs() < 1e-12);
        assert!(actual < vapor_pressure_pa);
    }

    #[test]
    fn mixture_threshold_is_bounded_by_population_extrema() {
        let result = evaluate_selective_cavitation_thresholds(&SelectiveCavitationInput {
            base_vapor_pressure_pa: 3_170.0,
            ambient_pressure_pa: 101_325.0,
            density_kg_m3: 1_000.0,
            populations: vec![healthy_population(), target_population()],
        })
        .expect("selective cavitation should evaluate");

        let min_threshold = result
            .population_thresholds
            .iter()
            .map(|entry| entry.effective_threshold_pressure_pa)
            .fold(f64::INFINITY, f64::min);
        let max_threshold = result
            .population_thresholds
            .iter()
            .map(|entry| entry.effective_threshold_pressure_pa)
            .fold(f64::NEG_INFINITY, f64::max);

        assert!(result.mixture_inception_threshold_pa >= min_threshold);
        assert!(result.mixture_inception_threshold_pa <= max_threshold);
    }

    proptest! {
        #[test]
        fn softer_larger_targets_raise_effective_threshold(
            healthy_sigma in 0.05_f64..0.09,
            target_sigma in 0.01_f64..0.04,
            healthy_radius in 3e-6_f64..6e-6,
            target_radius in 7e-6_f64..15e-6
        ) {
            let healthy = SelectiveCavitationPopulation {
                mechanical_state: CellMechanicalState {
                    membrane_stiffness_pa: 120_000.0,
                    interfacial_tension_n_m: healthy_sigma,
                    particle_radius_m: healthy_radius,
                    deformability_factor: 1.0,
                },
                ..healthy_population()
            };
            let target = SelectiveCavitationPopulation {
                mechanical_state: CellMechanicalState {
                    membrane_stiffness_pa: 20_000.0,
                    interfacial_tension_n_m: target_sigma,
                    particle_radius_m: target_radius,
                    deformability_factor: 1.2,
                },
                ..target_population()
            };

            let result = evaluate_selective_cavitation_thresholds(&SelectiveCavitationInput {
                base_vapor_pressure_pa: 3_170.0,
                ambient_pressure_pa: 101_325.0,
                density_kg_m3: 1_000.0,
                populations: vec![healthy, target],
            }).expect("selective cavitation should evaluate");

            prop_assert!(result.selectivity_margin_pa > 0.0);
        }
    }
}
