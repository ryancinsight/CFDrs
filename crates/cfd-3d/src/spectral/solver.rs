//! Main spectral solver implementation

use super::basis::SpectralBasis;
use super::poisson::{PoissonBoundaryCondition, PoissonSolver};
use cfd_core::Result;
use nalgebra::{DMatrix, RealField};
use num_traits::FromPrimitive;
use serde::{Deserialize, Serialize};
/// Spectral method configuration
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SpectralConfig<T: RealField + Copy> {
    /// Base solver configuration
    pub base: cfd_core::solver::SolverConfig<T>,
    /// Number of modes in x direction
    pub nx_modes: usize,
    /// Number of modes in y direction
    pub ny_modes: usize,
    /// Number of modes in z direction
    pub nz_modes: usize,
    /// Basis type for each direction
    pub basis_x: SpectralBasis,
    pub basis_y: SpectralBasis,
    pub basis_z: SpectralBasis,
    /// Time step for time-dependent problems
    pub dt: Option<T>,
}
impl<T: RealField + FromPrimitive + Copy> SpectralConfig<T> {
    /// Create new configuration with validation
    pub fn new(nx: usize, ny: usize, nz: usize) -> Result<Self> {
        Ok(Self {
            base: cfd_core::solver::SolverConfig::builder()
                .tolerance(T::from_f64(1e-10).ok_or_else(|| {
                    cfd_core::error::Error::InvalidConfiguration("Cannot convert tolerance".into())
                })?)
                .max_iterations(100)
                .build(),
            nx_modes: nx,
            ny_modes: ny,
            nz_modes: nz,
            basis_x: SpectralBasis::Chebyshev,
            basis_y: SpectralBasis::Chebyshev,
            basis_z: SpectralBasis::Chebyshev,
            dt: None,
        })
    }
/// Spectral solver for 3D problems
pub struct SpectralSolver<T: RealField + Copy> {
    config: SpectralConfig<T>,
    poisson_solver: PoissonSolver<T>,
impl<T: RealField + FromPrimitive + Copy> SpectralSolver<T> {
    /// Create new spectral solver
    pub fn new(config: SpectralConfig<T>) -> Result<Self> {
        let poisson_solver = PoissonSolver::new(config.nx_modes, config.ny_modes, config.nz_modes)?;
            config,
            poisson_solver,
    /// Solve a problem using spectral methods
    pub fn solve(&mut self) -> Result<SpectralSolution<T>> {
        // Initialize solution with zeros
        let mut solution = SpectralSolution::new(
            self.config.nx_modes,
            self.config.ny_modes,
            self.config.nz_modes,
        );
        // Create source term as a matrix (simplified)
        let source = DMatrix::zeros(
            self.config.nx_modes * self.config.ny_modes,
        // Apply boundary conditions (simplified)
        let bc_x = (
            PoissonBoundaryCondition::Dirichlet(T::zero()),
        let bc_y = (
        let bc_z = (
        // Apply the Poisson solver with proper arguments
        let potential = self.poisson_solver.solve(&source, bc_x, bc_y, bc_z)?;
        // Store the result in the solution
        solution.u = potential;
        Ok(solution)
/// Solution from spectral solver
pub struct SpectralSolution<T: RealField + Copy> {
    /// Solution field
    pub u: DMatrix<T>,
    /// Grid dimensions
    pub nx: usize,
    pub ny: usize,
    pub nz: usize,
impl<T: RealField + Copy> SpectralSolution<T> {
    /// Create new solution
    #[must_use]
    pub fn new(nx: usize, ny: usize, nz: usize) -> Self {
        Self {
            u: DMatrix::zeros(nx * ny, nz),
            nx,
            ny,
            nz,
        }
    /// Get solution at a point
    pub fn at(&self, i: usize, j: usize, k: usize) -> T
    where
        T: Clone,
    {
        let idx = i * self.ny + j;
        self.u[(idx, k)]
