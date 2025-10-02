//! Core traits for collision operators

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
