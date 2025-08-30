//! Core boundary condition types and fundamental classifications
//!
//! This module defines the basic boundary condition types used throughout
//! the CFD solver, following standard classifications from computational
//! fluid dynamics literature.
//!
//! Reference: Versteeg & Malalasekera (2007). An Introduction to Computational
//! Fluid Dynamics: The Finite Volume Method, Chapter 11.

use nalgebra::RealField;
use serde::{Deserialize, Serialize};

/// Fundamental physical classification of boundary conditions
///
/// Reference: LeVeque, R.J. (2002). Finite Volume Methods for Hyperbolic Problems, p. 213
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FundamentalBCType {
    /// Fixed value boundary condition (first type)
    Dirichlet,
    /// Fixed gradient boundary condition (second type)
    Neumann,
    /// Mixed boundary condition (third type)
    Robin,
    /// Other specialized boundary conditions
    Other,
}

/// Basic boundary condition types
///
/// These are the fundamental mathematical boundary conditions that can be
/// applied to partial differential equations.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub enum MathematicalBoundaryCondition<T: RealField + Copy> {
    /// Dirichlet: u = g on boundary
    Dirichlet { value: T },
    
    /// Neumann: ∂u/∂n = g on boundary
    Neumann { gradient: T },
    
    /// Robin: αu + β∂u/∂n = γ on boundary
    Robin { alpha: T, beta: T, gamma: T },
    
    /// Periodic: u(x) = u(x + L)
    Periodic { partner: String },
}

impl<T: RealField + Copy> MathematicalBoundaryCondition<T> {
    /// Returns the fundamental type classification
    pub const fn fundamental_type(&self) -> FundamentalBCType {
        match self {
            Self::Dirichlet { .. } => FundamentalBCType::Dirichlet,
            Self::Neumann { .. } => FundamentalBCType::Neumann,
            Self::Robin { .. } => FundamentalBCType::Robin,
            Self::Periodic { .. } => FundamentalBCType::Other,
        }
    }
    
    /// Check if this is a Dirichlet-type boundary
    pub const fn is_dirichlet(&self) -> bool {
        matches!(self, Self::Dirichlet { .. })
    }
    
    /// Check if this is a Neumann-type boundary
    pub const fn is_neumann(&self) -> bool {
        matches!(self, Self::Neumann { .. })
    }
}