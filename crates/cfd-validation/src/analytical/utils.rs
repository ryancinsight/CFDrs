//! Utility functions for analytical solutions

use crate::scalar;
use eunomia::FloatElement;
use eunomia::RealField;

/// Analytical utilities for common calculations
pub struct AnalyticalUtils;

impl AnalyticalUtils {
    /// Calculate Reynolds number
    pub fn reynolds_number<T: RealField + Copy>(
        density: T,
        velocity: T,
        length: T,
        viscosity: T,
    ) -> T {
        density * velocity * length / viscosity
    }

    /// Calculate Peclet number
    pub fn peclet_number<T: RealField + Copy>(velocity: T, length: T, diffusivity: T) -> T {
        velocity * length / diffusivity
    }

    /// Calculate Strouhal number
    pub fn strouhal_number<T: RealField + Copy>(frequency: T, length: T, velocity: T) -> T {
        frequency * length / velocity
    }

    /// Calculate Froude number
    pub fn froude_number<T: RealField + Copy + FloatElement>(
        velocity: T,
        length: T,
        gravity: T,
    ) -> T {
        velocity / scalar::sqrt(gravity * length)
    }

    /// Calculate Mach number
    pub fn mach_number<T: RealField + Copy>(velocity: T, sound_speed: T) -> T {
        velocity / sound_speed
    }

    /// Calculate pressure coefficient
    pub fn pressure_coefficient<T: RealField + Copy + FloatElement>(
        pressure: T,
        reference_pressure: T,
        density: T,
        reference_velocity: T,
    ) -> T {
        let half = scalar::from_f64::<T>(0.5);
        let dynamic_pressure = half * density * reference_velocity * reference_velocity;
        (pressure - reference_pressure) / dynamic_pressure
    }

    /// Calculate skin friction coefficient
    pub fn skin_friction_coefficient<T: RealField + Copy + FloatElement>(
        wall_shear_stress: T,
        density: T,
        reference_velocity: T,
    ) -> T {
        let half = scalar::from_f64::<T>(0.5);
        let dynamic_pressure = half * density * reference_velocity * reference_velocity;
        wall_shear_stress / dynamic_pressure
    }

    /// Calculate drag coefficient
    pub fn drag_coefficient<T: RealField + Copy + FloatElement>(
        drag_force: T,
        density: T,
        velocity: T,
        reference_area: T,
    ) -> T {
        let half = scalar::from_f64::<T>(0.5);
        let dynamic_pressure = half * density * velocity * velocity;
        drag_force / (dynamic_pressure * reference_area)
    }

    /// Calculate lift coefficient
    pub fn lift_coefficient<T: RealField + Copy + FloatElement>(
        lift_force: T,
        density: T,
        velocity: T,
        reference_area: T,
    ) -> T {
        let half = scalar::from_f64::<T>(0.5);
        let dynamic_pressure = half * density * velocity * velocity;
        lift_force / (dynamic_pressure * reference_area)
    }

    /// Check if flow is laminar based on Reynolds number
    #[allow(clippy::match_same_arms)] // Same critical Re for Sphere/Cylinder per fluid mechanics
    pub fn is_laminar<T: RealField + Copy + FloatElement>(
        reynolds: T,
        geometry: FlowGeometry,
    ) -> bool {
        let critical_re = match geometry {
            FlowGeometry::Pipe => scalar::from_f64(2300.0),
            FlowGeometry::FlatPlate => scalar::from_f64(500_000.0),
            FlowGeometry::Sphere | FlowGeometry::Cylinder => scalar::from_f64(200_000.0),
        };
        reynolds < critical_re
    }
}

/// Flow geometry type for Reynolds number interpretation
#[derive(Debug, Clone, Copy)]
pub enum FlowGeometry {
    /// Pipe flow
    Pipe,
    /// Flow over flat plate
    FlatPlate,
    /// Flow around sphere
    Sphere,
    /// Flow around cylinder
    Cylinder,
}
