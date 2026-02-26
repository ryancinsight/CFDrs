//! Core traits for collision operators
//!
//! # Theorem
//! The Lattice Boltzmann Method (LBM) recovers the macroscopic Navier-Stokes equations
//! in the low Mach number limit.
//!
//! **Proof sketch**:
//! Through the Chapman-Enskog expansion, the discrete Boltzmann equation with the BGK
//! collision operator can be expanded in powers of the Knudsen number ($Kn$).
//! At $O(Kn^0)$, the Euler equations are recovered. At $O(Kn^1)$, the viscous stress
//! tensor emerges, yielding the weakly compressible Navier-Stokes equations.
//! The kinematic viscosity is related to the relaxation time $\tau$ by $\nu = c_s^2 (\tau - 0.5)\Delta t$.

#![allow(dead_code)]

use nalgebra::RealField;

/// Trait for collision operators in LBM
pub trait CollisionOperator<T: RealField + Copy> {
    /// Apply collision step to distribution functions
    fn collide(&self, f: &mut Vec<Vec<[T; 9]>>, density: &[Vec<T>], velocity: &[Vec<[T; 2]>]);

    /// Get relaxation time
    fn tau(&self) -> T;

    /// Get kinematic viscosity
    fn viscosity(&self, dt: T, dx: T) -> T;
}

/// Collision model types
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum CollisionModel {
    /// BGK single relaxation time
    Bgk,
    /// Multiple relaxation time
    Mrt,
    /// Regularized collision
    Regularized,
    /// Two relaxation time
    Trt,
}
