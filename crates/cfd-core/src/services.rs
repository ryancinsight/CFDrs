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
    pub fn reynolds_number<T: RealField + num_traits::Float>(
        fluid: &Fluid<T>,
        velocity: T,
        characteristic_length: T,
    ) -> T {
        velocity * characteristic_length / fluid.kinematic_viscosity()
    }

    /// Calculate Prandtl number if thermal properties are available
    pub fn prandtl_number<T: RealField + num_traits::Float>(fluid: &Fluid<T>) -> Option<T> {
        match (fluid.specific_heat.clone(), fluid.thermal_conductivity.clone()) {
            (Some(cp), Some(k)) => Some(fluid.characteristic_viscosity() * cp / k),
            _ => None,
        }
    }

    /// Determine flow regime based on Reynolds number
    pub fn flow_regime<T: RealField + FromPrimitive>(reynolds: T) -> FlowRegime {
        let re_2300 = T::from_f64(2300.0).ok_or_else(|| crate::error::Error::Numerical(crate::error::NumericalErrorKind::InvalidFpOperation))?;
        let re_4000 = T::from_f64(4000.0).ok_or_else(|| crate::error::Error::Numerical(crate::error::NumericalErrorKind::InvalidFpOperation))?;

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
        let re_2300 = T::from_f64(2300.0).ok_or_else(|| crate::error::Error::Numerical(crate::error::NumericalErrorKind::InvalidFpOperation))?;
        let sixty_four = T::from_f64(64.0).ok_or_else(|| crate::error::Error::Numerical(crate::error::NumericalErrorKind::InvalidFpOperation))?;

        if reynolds < re_2300 {
            // Laminar flow: f = 64/Re
            Ok(sixty_four / reynolds)
        } else {
            // Turbulent flow: use Colebrook-White or smooth pipe approximation
            match roughness {
                Some(eps) => {
                    // Colebrook-White equation
                    let relative_roughness = eps / diameter;
                    Self::colebrook_white_friction_factor(reynolds, relative_roughness)
                }
                None => {
                    // Smooth pipe: Blasius equation for Re < 100,000
                    let re_turbulent = T::from_f64(100_000.0).ok_or_else(|| crate::error::Error::Numerical(crate::error::NumericalErrorKind::InvalidFpOperation))?;
                    if reynolds < re_turbulent {
                        let coeff = T::from_f64(0.316).ok_or_else(|| crate::error::Error::Numerical(crate::error::NumericalErrorKind::InvalidFpOperation))?;
                        let exp = T::from_f64(0.25).ok_or_else(|| crate::error::Error::Numerical(crate::error::NumericalErrorKind::InvalidFpOperation))?;
                        Ok(coeff / reynolds.powf(exp))
                    } else {
                        // Use Prandtl-Karman equation
                        Self::prandtl_karman_friction_factor(reynolds)
                    }
                }
            }
        }
    }

    /// Colebrook-White friction factor (iterative solution)
    fn colebrook_white_friction_factor<T: RealField + FromPrimitive + Copy>(
        reynolds: T,
        relative_roughness: T,
    ) -> Result<T> {
        let mut f = T::from_f64(0.02).ok_or_else(|| crate::error::Error::Numerical(crate::error::NumericalErrorKind::InvalidFpOperation))?; // Initial guess
        let tolerance = T::from_f64(1e-6).ok_or_else(|| crate::error::Error::Numerical(crate::error::NumericalErrorKind::InvalidFpOperation))?;
        let max_iterations = 50;

        for _ in 0..max_iterations {
            let term1 = relative_roughness / T::from_f64(3.7).ok_or_else(|| crate::error::Error::Numerical(crate::error::NumericalErrorKind::InvalidFpOperation))?;
            let term2 = T::from_f64(2.51).ok_or_else(|| crate::error::Error::Numerical(crate::error::NumericalErrorKind::InvalidFpOperation))? / (reynolds * f.sqrt());
            let f_new = T::from_f64(0.25).ok_or_else(|| crate::error::Error::Numerical(crate::error::NumericalErrorKind::InvalidFpOperation))? / 
                (T::from_f64(10.0).ok_or_else(|| crate::error::Error::Numerical(crate::error::NumericalErrorKind::InvalidFpOperation))?.ln() * (term1 + term2)).powi(2);
            
            if (f_new - f).abs() < tolerance {
                return Ok(f_new);
            }
            f = f_new;
        }

        Err(Error::Convergence(crate::error::ConvergenceErrorKind::MaxIterationsExceeded { max: 100 }
        ))
    }

    /// Prandtl-Karman friction factor for smooth pipes
    fn prandtl_karman_friction_factor<T: RealField + FromPrimitive + Copy>(reynolds: T) -> Result<T> {
        let log_re = reynolds.ln() / T::from_f64(10.0).ok_or_else(|| crate::error::Error::Numerical(crate::error::NumericalErrorKind::InvalidFpOperation))?.ln();
        let coeff = T::from_f64(2.0).ok_or_else(|| crate::error::Error::Numerical(crate::error::NumericalErrorKind::InvalidFpOperation))? * log_re - T::from_f64(0.8).ok_or_else(|| crate::error::Error::Numerical(crate::error::NumericalErrorKind::InvalidFpOperation))?;
        Ok(T::one() / (coeff * coeff))
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
    pub fn assess_quality<T: RealField + FromPrimitive>(
        aspect_ratio_stats: &QualityStatistics<T>,
        skewness_stats: &QualityStatistics<T>,
        orthogonality_stats: &QualityStatistics<T>,
    ) -> QualityAssessment {
        let aspect_ratio_quality = Self::assess_aspect_ratio(aspect_ratio_stats);
        let skewness_quality = Self::assess_skewness(skewness_stats);
        let orthogonality_quality = Self::assess_orthogonality(orthogonality_stats);

        let overall_quality = [aspect_ratio_quality, skewness_quality, orthogonality_quality]
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

    fn assess_aspect_ratio<T: RealField + FromPrimitive>(stats: &QualityStatistics<T>) -> QualityLevel {
        let threshold_level4 = T::from_f64(2.0).ok_or_else(|| crate::error::Error::Numerical(crate::error::NumericalErrorKind::InvalidFpOperation))?;
        let threshold_level3 = T::from_f64(5.0).ok_or_else(|| crate::error::Error::Numerical(crate::error::NumericalErrorKind::InvalidFpOperation))?;
        let threshold_level2 = T::from_f64(10.0).ok_or_else(|| crate::error::Error::Numerical(crate::error::NumericalErrorKind::InvalidFpOperation))?;

        if stats.max < threshold_level4 {
            QualityLevel::Level4
        } else if stats.max < threshold_level3 {
            QualityLevel::Level3
        } else if stats.max < threshold_level2 {
            QualityLevel::Level2
        } else {
            QualityLevel::Level1
        }
    }

    fn assess_skewness<T: RealField + FromPrimitive>(stats: &QualityStatistics<T>) -> QualityLevel {
        let threshold_level4 = T::from_f64(0.25).ok_or_else(|| crate::error::Error::Numerical(crate::error::NumericalErrorKind::InvalidFpOperation))?;
        let threshold_level3 = T::from_f64(0.5).ok_or_else(|| crate::error::Error::Numerical(crate::error::NumericalErrorKind::InvalidFpOperation))?;
        let threshold_level2 = T::from_f64(0.8).ok_or_else(|| crate::error::Error::Numerical(crate::error::NumericalErrorKind::InvalidFpOperation))?;

        if stats.max < threshold_level4 {
            QualityLevel::Level4
        } else if stats.max < threshold_level3 {
            QualityLevel::Level3
        } else if stats.max < threshold_level2 {
            QualityLevel::Level2
        } else {
            QualityLevel::Level1
        }
    }

    fn assess_orthogonality<T: RealField + FromPrimitive>(stats: &QualityStatistics<T>) -> QualityLevel {
        let threshold_level4 = T::from_f64(0.95).ok_or_else(|| crate::error::Error::Numerical(crate::error::NumericalErrorKind::InvalidFpOperation))?;
        let threshold_level3 = T::from_f64(0.85).ok_or_else(|| crate::error::Error::Numerical(crate::error::NumericalErrorKind::InvalidFpOperation))?;
        let threshold_level2 = T::from_f64(0.7).ok_or_else(|| crate::error::Error::Numerical(crate::error::NumericalErrorKind::InvalidFpOperation))?;

        if stats.min > threshold_level4 {
            QualityLevel::Level4
        } else if stats.min > threshold_level3 {
            QualityLevel::Level3
        } else if stats.min > threshold_level2 {
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
            recommendations.push("Consider refining mesh in regions with high aspect ratio cells".to_string());
        }

        if skewness == QualityLevel::Level1 {
            recommendations.push("Improve mesh quality by reducing cell skewness".to_string());
        }

        if orthogonality == QualityLevel::Level1 {
            recommendations.push("Enhance mesh orthogonality for better numerical accuracy".to_string());
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
pub struct QualityStatistics<T: RealField> {
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
