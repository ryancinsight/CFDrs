//! Utility functions for analytical solutions

use nalgebra::RealField;
use num_traits::FromPrimitive;

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
    pub fn froude_number<T: RealField + Copy + FromPrimitive>(
        velocity: T,
        length: T,
        gravity: T,
    ) -> T {
        velocity / (gravity * length).sqrt()
    }

    /// Calculate Mach number
    pub fn mach_number<T: RealField + Copy>(velocity: T, sound_speed: T) -> T {
        velocity / sound_speed
    }

    /// Calculate pressure coefficient
    pub fn pressure_coefficient<T: RealField + Copy + FromPrimitive>(
        pressure: T,
        reference_pressure: T,
        density: T,
        reference_velocity: T,
    ) -> T {
        let half = T::from_f64(0.5).unwrap_or(T::one() / (T::one() + T::one()));
        let dynamic_pressure = half * density * reference_velocity * reference_velocity;
        (pressure - reference_pressure) / dynamic_pressure
    }

    /// Calculate skin friction coefficient
    pub fn skin_friction_coefficient<T: RealField + Copy + FromPrimitive>(
        wall_shear_stress: T,
        density: T,
        reference_velocity: T,
    ) -> T {
        let half = T::from_f64(0.5).unwrap_or(T::one() / (T::one() + T::one()));
        let dynamic_pressure = half * density * reference_velocity * reference_velocity;
        wall_shear_stress / dynamic_pressure
    }

    /// Calculate drag coefficient
    pub fn drag_coefficient<T: RealField + Copy + FromPrimitive>(
        drag_force: T,
        density: T,
        velocity: T,
        reference_area: T,
    ) -> T {
        let half = T::from_f64(0.5).unwrap_or(T::one() / (T::one() + T::one()));
        let dynamic_pressure = half * density * velocity * velocity;
        drag_force / (dynamic_pressure * reference_area)
    }

    /// Calculate lift coefficient
    pub fn lift_coefficient<T: RealField + Copy + FromPrimitive>(
        lift_force: T,
        density: T,
        velocity: T,
        reference_area: T,
    ) -> T {
        let half = T::from_f64(0.5).unwrap_or(T::one() / (T::one() + T::one()));
        let dynamic_pressure = half * density * velocity * velocity;
        lift_force / (dynamic_pressure * reference_area)
    }

    /// Check if flow is laminar based on Reynolds number
    pub fn is_laminar<T: RealField + Copy + FromPrimitive>(
        reynolds: T,
        geometry: FlowGeometry,
    ) -> bool {
        let critical_re = match geometry {
            FlowGeometry::Pipe => {
                T::from_f64(2300.0).unwrap_or(T::from_f64(2000.0).unwrap_or(T::one()))
            }
            FlowGeometry::FlatPlate => {
                T::from_f64(500000.0).unwrap_or(T::from_f64(100000.0).unwrap_or(T::one()))
            }
            FlowGeometry::Sphere => {
                T::from_f64(200000.0).unwrap_or(T::from_f64(100000.0).unwrap_or(T::one()))
            }
            FlowGeometry::Cylinder => {
                T::from_f64(200000.0).unwrap_or(T::from_f64(100000.0).unwrap_or(T::one()))
            }
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
