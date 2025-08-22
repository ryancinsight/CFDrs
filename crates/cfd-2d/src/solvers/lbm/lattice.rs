//! Lattice models for LBM simulations.
//!
//! This module defines lattice structures and their properties,
//! including velocity sets, weights, and related constants.

use nalgebra::RealField;

/// Trait defining a lattice model for LBM
pub trait LatticeModel {
    /// Number of discrete velocities
    const Q: usize;
    
    /// Get lattice velocities
    fn velocities() -> &'static [(i32, i32)];
    
    /// Get lattice weights
    fn weights() -> &'static [f64];
    
    /// Get opposite direction index
    fn opposite(direction: usize) -> usize;
}

/// D2Q9 lattice model (2D, 9 velocities)
pub struct D2Q9;

impl D2Q9 {
    /// Number of velocity directions
    pub const Q: usize = 9;

    /// Lattice velocities (normalized by lattice spacing)
    pub const VELOCITIES: [(i32, i32); 9] = [
        (0, 0),   // 0: rest
        (1, 0),   // 1: east
        (0, 1),   // 2: north
        (-1, 0),  // 3: west
        (0, -1),  // 4: south
        (1, 1),   // 5: northeast
        (-1, 1),  // 6: northwest
        (-1, -1), // 7: southwest
        (1, -1),  // 8: southeast
    ];

    /// Lattice weights
    pub const WEIGHTS: [f64; 9] = [
        4.0/9.0,  // 0: rest
        1.0/9.0,  // 1-4: cardinal directions
        1.0/9.0,
        1.0/9.0,
        1.0/9.0,
        1.0/36.0, // 5-8: diagonal directions
        1.0/36.0,
        1.0/36.0,
        1.0/36.0,
    ];

    /// Opposite direction indices
    pub const OPPOSITE: [usize; 9] = [0, 3, 4, 1, 2, 7, 8, 5, 6];
}

impl LatticeModel for D2Q9 {
    const Q: usize = 9;
    
    fn velocities() -> &'static [(i32, i32)] {
        &Self::VELOCITIES
    }
    
    fn weights() -> &'static [f64] {
        &Self::WEIGHTS
    }
    
    fn opposite(direction: usize) -> usize {
        Self::OPPOSITE[direction]
    }
}

/// Compute equilibrium distribution function
pub fn equilibrium<T: RealField + Copy>(
    density: T,
    velocity: &[T; 2],
    direction: usize,
    weight: T,
    lattice_velocity: &(i32, i32),
) -> T {
    let cx = T::from_i32(lattice_velocity.0).unwrap_or_else(T::zero);
    let cy = T::from_i32(lattice_velocity.1).unwrap_or_else(T::zero);
    
    let u_sq = velocity[0] * velocity[0] + velocity[1] * velocity[1];
    let cu = cx * velocity[0] + cy * velocity[1];
    
    let three = T::from_f64(3.0).unwrap_or_else(T::zero);
    let nine_half = T::from_f64(4.5).unwrap_or_else(T::zero);
    let three_half = T::from_f64(1.5).unwrap_or_else(T::zero);
    
    weight * density * (T::one() + three * cu + nine_half * cu * cu - three_half * u_sq)
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_d2q9_properties() {
        assert_eq!(D2Q9::Q, 9);
        assert_eq!(D2Q9::VELOCITIES.len(), 9);
        assert_eq!(D2Q9::WEIGHTS.len(), 9);
        assert_eq!(D2Q9::OPPOSITE.len(), 9);
        
        // Test weight sum equals 1
        let weight_sum: f64 = D2Q9::WEIGHTS.iter().sum();
        assert!((weight_sum - 1.0).abs() < 1e-10);
        
        // Test opposite directions
        for i in 0..9 {
            let opp = D2Q9::opposite(i);
            assert_eq!(D2Q9::opposite(opp), i);
        }
    }
}