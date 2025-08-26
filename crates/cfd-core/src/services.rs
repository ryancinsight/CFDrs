//! Domain services for CFD computations.
//!
//! This module provides domain services that encapsulate complex business logic
//! and coordinate between different domain entities following DDD principles.

use crate::error::{Error, Result};
use crate::fluid::Fluid;
use nalgebra::RealField;
use num_traits::FromPrimitive;

/// Service for fluid dynamics calculations
pub struct FluidDynamicsService;

impl FluidDynamicsService {
    /// Calculate Reynolds number for a given flow configuration
    pub fn reynolds_number<T: RealField + Copy + num_traits::Float>(
        fluid: &Fluid<T>,
        velocity: T,
        characteristic_length: T,
    ) -> T {
        velocity * characteristic_length / fluid.kinematic_viscosity()
    }

    /// Calculate Prandtl number if thermal properties are available
    pub fn prandtl_number<T: RealField + Copy + num_traits::Float>(fluid: &Fluid<T>) -> Option<T> {
        match (fluid.specific_heat, fluid.thermal_conductivity) {
            (Some(cp), Some(k)) => Some(fluid.characteristic_viscosity() * cp / k),
            _ => None,
        }
    }

    /// Determine flow regime based on Reynolds number
    pub fn flow_regime<T: RealField + FromPrimitive + Copy>(reynolds: T) -> FlowRegime {
        let re_2300 = T::from_f64(crate::constants::dimensionless::reynolds::PIPE_CRITICAL_LOWER)
            .unwrap_or_else(|| T::one());
        let re_4000 = T::from_f64(crate::constants::dimensionless::reynolds::PIPE_CRITICAL_UPPER)
            .unwrap_or_else(|| T::one());

        if reynolds < re_2300 {
            FlowRegime::Laminar
        } else if reynolds > re_4000 {
            FlowRegime::Turbulent
        } else {
            FlowRegime::Transitional
        }
    }

    /// Calculate pressure drop for pipe flow
    pub fn pipe_pressure_drop<T: RealField + FromPrimitive + Copy + num_traits::Float>(
        fluid: &Fluid<T>,
        velocity: T,
        length: T,
        diameter: T,
        roughness: Option<T>,
    ) -> Result<T> {
        let reynolds = Self::reynolds_number(fluid, velocity, diameter);
        let friction_factor = Self::friction_factor(reynolds, diameter, roughness)?;

        let two = T::one() + T::one();

        Ok(friction_factor * length * fluid.density * velocity * velocity / (two * diameter))
    }

    /// Calculate friction factor using appropriate correlation
    fn friction_factor<T: RealField + FromPrimitive + Copy>(
        reynolds: T,
        diameter: T,
        roughness: Option<T>,
    ) -> Result<T> {
        let re_2300 = T::from_f64(crate::constants::dimensionless::reynolds::PIPE_CRITICAL_LOWER)
            .unwrap_or_else(|| T::one());
        let sixty_four = T::from_f64(64.0).unwrap_or_else(|| T::one());

        if reynolds < re_2300 {
            // Laminar flow: f = 64/Re
            Ok(sixty_four / reynolds)
        } else {
            // Turbulent flow: use Colebrook-White or explicit correlations
            if let Some(eps) = roughness {
                // Colebrook-White equation
                let relative_roughness = eps / diameter;
                Self::colebrook_white_friction_factor(reynolds, relative_roughness)
            } else {
                // Smooth pipe: Blasius (low Re) or Haaland (general explicit)
                let blasius_max =
                    T::from_f64(crate::constants::physics::hydraulics::BLASIUS_MAX_RE)
                        .unwrap_or_else(|| T::one());
                if reynolds < blasius_max {
                    let coeff =
                        T::from_f64(crate::constants::physics::hydraulics::BLASIUS_COEFFICIENT)
                            .unwrap_or_else(|| T::one());
                    let exp = T::from_f64(crate::constants::physics::hydraulics::BLASIUS_EXPONENT)
                        .unwrap_or_else(|| T::one());
                    Ok(coeff / reynolds.powf(exp))
                } else {
                    Self::haaland_friction_factor(reynolds, T::zero())
                }
            }
        }
    }

    /// Colebrook-White friction factor (iterative solution)
    fn colebrook_white_friction_factor<T: RealField + FromPrimitive + Copy>(
        reynolds: T,
        relative_roughness: T,
    ) -> Result<T> {
        let mut f = T::from_f64(0.02).unwrap_or_else(|| T::one());
        let tolerance = T::from_f64(1e-6).unwrap_or_else(|| T::one());
        let max_iterations = 50;

        for _ in 0..max_iterations {
            // 1/sqrt(f) = -2 log10( (ε/D)/3.7 + 2.51/(Re sqrt(f)) )
            // Use log10(x) = ln(x) / ln(10)
            let divisor =
                T::from_f64(crate::constants::physics::hydraulics::COLEBROOK_ROUGHNESS_DIVISOR)
                    .unwrap_or_else(|| T::one());
            let num =
                T::from_f64(crate::constants::physics::hydraulics::COLEBROOK_REYNOLDS_NUMERATOR)
                    .unwrap_or_else(|| T::one());
            let ln10 = T::from_f64(10.0).unwrap_or_else(|| T::one()).ln();

            let term1 = relative_roughness / divisor;
            let term2 = num / (reynolds * f.sqrt());
            let inside = term1 + term2;
            // Guard against invalid log argument
            if inside <= T::zero() {
                return Err(Error::Numerical(
                    crate::error::NumericalErrorKind::InvalidValue {
                        value: "Colebrook-White log argument".to_string(),
                    },
                ));
            }
            let inv_sqrt_f = -(T::from_f64(2.0).unwrap_or_else(|| T::one())) * (inside.ln() / ln10);
            let f_new = T::one() / (inv_sqrt_f * inv_sqrt_f);

            if (f_new - f).abs() < tolerance {
                return Ok(f_new);
            }
            f = f_new;
        }

        Err(Error::Convergence(
            crate::error::ConvergenceErrorKind::MaxIterationsExceeded { max: 100 },
        ))
    }

    /// Haaland explicit friction factor correlation
    fn haaland_friction_factor<T: RealField + FromPrimitive + Copy>(
        reynolds: T,
        relative_roughness: T,
    ) -> Result<T> {
        let a = T::from_f64(crate::constants::physics::hydraulics::HAALAND_LOG_COEFFICIENT)
            .unwrap_or_else(|| T::one());
        let rr_exp = T::from_f64(crate::constants::physics::hydraulics::HAALAND_ROUGHNESS_EXPONENT)
            .unwrap_or_else(|| T::one());
        let divisor =
            T::from_f64(crate::constants::physics::hydraulics::COLEBROOK_ROUGHNESS_DIVISOR)
                .unwrap_or_else(|| T::one());
        let ln10 = T::from_f64(10.0).unwrap_or_else(|| T::one()).ln();

        let term_rr = (relative_roughness / divisor).powf(rr_exp);
        let term_re = T::from_f64(6.9).unwrap_or_else(|| T::one()) / reynolds;
        let inside = term_rr + term_re;
        if inside <= T::zero() {
            return Err(Error::Numerical(
                crate::error::NumericalErrorKind::InvalidValue {
                    value: "Haaland log argument".to_string(),
                },
            ));
        }
        let inv_sqrt_f = a * (inside.ln() / ln10);
        Ok(T::one() / (inv_sqrt_f * inv_sqrt_f))
    }
}

/// Flow regime classification
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FlowRegime {
    /// Laminar flow (Re < 2300)
    Laminar,
    /// Transitional flow (2300 < Re < 4000)
    Transitional,
    /// Turbulent flow (Re > 4000)
    Turbulent,
}

/// Service for mesh quality assessment
pub struct MeshQualityService;

impl MeshQualityService {
    /// Assess overall mesh quality and provide recommendations
    pub fn assess_quality<T: RealField + FromPrimitive + Copy>(
        aspect_ratio_stats: &QualityStatistics<T>,
        skewness_stats: &QualityStatistics<T>,
        orthogonality_stats: &QualityStatistics<T>,
    ) -> QualityAssessment {
        let aspect_ratio_quality = Self::assess_aspect_ratio(aspect_ratio_stats);
        let skewness_quality = Self::assess_skewness(skewness_stats);
        let orthogonality_quality = Self::assess_orthogonality(orthogonality_stats);

        let overall_quality = [
            aspect_ratio_quality,
            skewness_quality,
            orthogonality_quality,
        ]
        .iter()
        .min()
        .copied()
        .unwrap_or(QualityLevel::Level1);

        QualityAssessment {
            overall_quality,
            aspect_ratio_quality,
            skewness_quality,
            orthogonality_quality,
            recommendations: Self::generate_recommendations(
                aspect_ratio_quality,
                skewness_quality,
                orthogonality_quality,
            ),
        }
    }

    fn assess_aspect_ratio<T: RealField + FromPrimitive + Copy>(
        stats: &QualityStatistics<T>,
    ) -> QualityLevel {
        let threshprevious_level4 = T::from_f64(2.0).unwrap_or_else(|| T::one());
        let threshprevious_level3 = T::from_f64(5.0).unwrap_or_else(|| T::one());
        let threshprevious_level2 = T::from_f64(10.0).unwrap_or_else(|| T::one());

        if stats.max < threshprevious_level4 {
            QualityLevel::Level4
        } else if stats.max < threshprevious_level3 {
            QualityLevel::Level3
        } else if stats.max < threshprevious_level2 {
            QualityLevel::Level2
        } else {
            QualityLevel::Level1
        }
    }

    fn assess_skewness<T: RealField + FromPrimitive + Copy>(
        stats: &QualityStatistics<T>,
    ) -> QualityLevel {
        let threshprevious_level4 = T::from_f64(0.25).unwrap_or_else(|| T::one());
        let threshprevious_level3 = T::from_f64(0.5).unwrap_or_else(|| T::one());
        let threshprevious_level2 = T::from_f64(0.8).unwrap_or_else(|| T::one());

        if stats.max < threshprevious_level4 {
            QualityLevel::Level4
        } else if stats.max < threshprevious_level3 {
            QualityLevel::Level3
        } else if stats.max < threshprevious_level2 {
            QualityLevel::Level2
        } else {
            QualityLevel::Level1
        }
    }

    fn assess_orthogonality<T: RealField + FromPrimitive + Copy>(
        stats: &QualityStatistics<T>,
    ) -> QualityLevel {
        let threshprevious_level4 = T::from_f64(0.95).unwrap_or_else(|| T::one());
        let threshprevious_level3 = T::from_f64(0.85).unwrap_or_else(|| T::one());
        let threshprevious_level2 = T::from_f64(0.7).unwrap_or_else(|| T::one());

        if stats.min > threshprevious_level4 {
            QualityLevel::Level4
        } else if stats.min > threshprevious_level3 {
            QualityLevel::Level3
        } else if stats.min > threshprevious_level2 {
            QualityLevel::Level2
        } else {
            QualityLevel::Level1
        }
    }

    fn generate_recommendations(
        aspect_ratio: QualityLevel,
        skewness: QualityLevel,
        orthogonality: QualityLevel,
    ) -> Vec<String> {
        let mut recommendations = Vec::new();

        if aspect_ratio == QualityLevel::Level1 {
            recommendations
                .push("Consider refining mesh in regions with high aspect ratio cells".to_string());
        }

        if skewness == QualityLevel::Level1 {
            recommendations.push("Improve mesh quality by reducing cell skewness".to_string());
        }

        if orthogonality == QualityLevel::Level1 {
            recommendations
                .push("Enhance mesh orthogonality for better numerical accuracy".to_string());
        }

        if recommendations.is_empty() {
            recommendations.push("Mesh quality is acceptable for CFD simulation".to_string());
        }

        recommendations
    }
}

/// Quality level enumeration
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub enum QualityLevel {
    /// Level 1 - may cause numerical issues
    Level1,
    /// Level 2 - usable for basic simulations
    Level2,
    /// Level 3 - suitable for most applications
    Level3,
    /// Level 4 - suitable for high-accuracy simulations
    Level4,
}

/// Quality statistics structure
#[derive(Debug, Clone)]
pub struct QualityStatistics<T: RealField + Copy> {
    /// Minimum value
    pub min: T,
    /// Maximum value
    pub max: T,
    /// Mean value
    pub mean: T,
    /// Standard deviation
    pub std_dev: T,
}

/// Overall quality assessment
#[derive(Debug, Clone)]
pub struct QualityAssessment {
    /// Overall mesh quality
    pub overall_quality: QualityLevel,
    /// Aspect ratio quality
    pub aspect_ratio_quality: QualityLevel,
    /// Skewness quality
    pub skewness_quality: QualityLevel,
    /// Orthogonality quality
    pub orthogonality_quality: QualityLevel,
    /// Improvement recommendations
    pub recommendations: Vec<String>,
}
