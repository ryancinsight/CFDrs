//! SUPG/PSPG stabilization for FEM
//!
//! References:
//! - Brooks, A.N. and Hughes, T.J.R. (1982). "Streamline upwind/Petrov-Galerkin formulations
//!   for convection dominated flows with particular emphasis on the incompressible Navier-Stokes equations"
//! - Tezduyar, T.E. (1991). "Stabilized finite element formulations for incompressible flow computations"

use nalgebra::{RealField, Vector3};
use num_traits::FromPrimitive;

// Named constants for stabilization
const TWO: f64 = 2.0;
const FOUR: f64 = 4.0;

/// SUPG/PSPG stabilization parameters
pub struct StabilizationParameters<T: RealField + Copy> {
    /// Element size
    h: T,
    /// Kinematic viscosity
    nu: T,
    /// Velocity magnitude
    u_mag: T,
    /// Time step (for transient problems)
    dt: Option<T>,
}

impl<T: RealField + FromPrimitive + Copy> StabilizationParameters<T> {
    /// Create new stabilization parameters
    pub fn new(h: T, nu: T, velocity: Vector3<T>, dt: Option<T>) -> Self {
        let u_mag = velocity.norm();
        Self { h, nu, u_mag, dt }
    }

    /// Calculate SUPG stabilization parameter tau
    ///
    /// Based on Tezduyar (1991) formulation:
    /// τ = [(2/Δt)² + (2U/h)² + (4ν/h²)²]^(-1/2)
    ///
    /// For steady-state: τ = [(2U/h)² + (4ν/h²)²]^(-1/2)
    pub fn tau_supg(&self) -> T {
        let two = T::from_f64(TWO).unwrap_or_else(T::zero);
        let four = T::from_f64(FOUR).unwrap_or_else(T::zero);

        // Advection term: (2U/h)²
        let advection_term = if self.u_mag > T::zero() {
            let term = (two * self.u_mag) / self.h;
            term * term
        } else {
            T::zero()
        };

        // Diffusion term: (4ν/h²)²
        let diffusion_term = {
            let term = (four * self.nu) / (self.h * self.h);
            term * term
        };

        // Time term: (2/Δt)² (only for transient)
        let time_term = if let Some(dt) = self.dt {
            let term = two / dt;
            term * term
        } else {
            T::zero()
        };

        // Combined tau
        let sum = time_term + advection_term + diffusion_term;
        if sum > T::zero() {
            T::one() / sum.sqrt()
        } else {
            T::zero()
        }
    }

    /// Calculate PSPG stabilization parameter tau
    ///
    /// For pressure stabilization, uses same formulation as SUPG
    pub fn tau_pspg(&self) -> T {
        self.tau_supg()
    }

    /// Calculate element Peclet number
    /// Pe = U*h/(2ν)
    pub fn peclet_number(&self) -> T {
        if self.nu > T::zero() {
            (self.u_mag * self.h) / (T::from_f64(TWO).unwrap_or_else(T::zero) * self.nu)
        } else {
            T::zero()
        }
    }

    /// Calculate element Reynolds number
    /// Re = U*h/ν
    pub fn element_reynolds(&self) -> T {
        if self.nu > T::zero() {
            (self.u_mag * self.h) / self.nu
        } else {
            T::zero()
        }
    }

    /// Get optimal stabilization based on flow regime
    pub fn optimal_tau(&self) -> T {
        let pe = self.peclet_number();
        let tau_supg = self.tau_supg();

        // Adjust based on Peclet number
        if pe < T::one() {
            // Diffusion-dominated: reduce stabilization
            tau_supg * pe
        } else if pe > T::from_f64(100.0).unwrap_or_else(T::zero) {
            // Highly advection-dominated: use full stabilization
            tau_supg
        } else {
            // Intermediate regime: smooth transition
            let factor = (T::one() - T::one() / pe).max(T::zero());
            tau_supg * factor
        }
    }
}

/// Calculate element size for different element types
pub fn calculate_element_size<T: RealField + FromPrimitive + Copy>(
    vertices: &[Vector3<T>],
    velocity_direction: &Vector3<T>,
) -> T {
    // For tetrahedral elements (4 vertices)
    if vertices.len() == 4 {
        calculate_tetrahedral_size(vertices, velocity_direction)
    } else if vertices.len() == 8 {
        // For hexahedral elements
        calculate_hexahedral_size(vertices, velocity_direction)
    } else {
        // Default: use minimum edge length
        calculate_min_edge_length(vertices)
    }
}

/// Calculate size for tetrahedral element
fn calculate_tetrahedral_size<T: RealField + FromPrimitive + Copy>(
    vertices: &[Vector3<T>],
    velocity_direction: &Vector3<T>,
) -> T {
    if velocity_direction.norm() > T::zero() {
        // Directional element size in flow direction
        let dir = velocity_direction.normalize();

        // Project element onto flow direction
        let mut min_proj = T::max_value().unwrap_or_else(T::one);
        let mut max_proj = T::min_value().unwrap_or_else(T::zero);

        for vertex in vertices {
            let proj = vertex.dot(&dir);
            min_proj = min_proj.min(proj);
            max_proj = max_proj.max(proj);
        }

        (max_proj - min_proj).abs()
    } else {
        // No flow: use characteristic length
        calculate_min_edge_length(vertices)
    }
}

/// Calculate size for hexahedral element
fn calculate_hexahedral_size<T: RealField + FromPrimitive + Copy>(
    vertices: &[Vector3<T>],
    velocity_direction: &Vector3<T>,
) -> T {
    // Similar to tetrahedral but for 8 vertices
    calculate_tetrahedral_size(vertices, velocity_direction)
}

/// Calculate minimum edge length of element
fn calculate_min_edge_length<T: RealField + Copy>(vertices: &[Vector3<T>]) -> T {
    let mut min_length = T::max_value().unwrap_or_else(T::one);

    for i in 0..vertices.len() {
        for j in i + 1..vertices.len() {
            let edge_length = (vertices[i] - vertices[j]).norm();
            min_length = min_length.min(edge_length);
        }
    }

    min_length
}
