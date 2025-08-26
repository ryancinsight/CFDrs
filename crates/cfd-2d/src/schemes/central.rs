//! Central difference schemes

use super::{Grid2D, SpatialDiscretization};
use nalgebra::RealField;
use num_traits::FromPrimitive;
/// Second-order central difference scheme
pub struct CentralDifference<T: RealField + Copy> {
    _phantom: std::marker::PhantomData<T>,
}
impl<T: RealField + Copy> Default for CentralDifference<T> {
    fn default() -> Self {
        Self::new()
    }
impl<T: RealField + Copy> CentralDifference<T> {
    /// Create new central difference scheme
    #[must_use]
    pub fn new() -> Self {
        Self {
            _phantom: std::marker::PhantomData,
        }
impl<T: RealField + Copy + FromPrimitive + Copy> SpatialDiscretization<T> for CentralDifference<T> {
    fn compute_derivative(&self, grid: &Grid2D<T>, i: usize, j: usize) -> T {
        let two = T::from_f64(2.0).unwrap_or_else(T::zero);
        (grid.data[(i + 1, j)] - grid.data[(i - 1, j)]) / (two * grid.dx)
    fn order(&self) -> usize {
        2
    fn is_conservative(&self) -> bool {
        true
/// Fourth-order central difference scheme
pub struct FourthOrderCentral<T: RealField + Copy> {
impl<T: RealField + Copy> Default for FourthOrderCentral<T> {
impl<T: RealField + Copy> FourthOrderCentral<T> {
    /// Create new fourth-order central scheme
impl<T: RealField + Copy + FromPrimitive + Copy> SpatialDiscretization<T>
    for FourthOrderCentral<T>
{
        let eight = T::from_f64(8.0).unwrap_or_else(T::zero);
        let twelve = T::from_f64(12.0).unwrap_or_else(T::zero);
        (-grid.data[(i + 2, j)] + eight * grid.data[(i + 1, j)] - eight * grid.data[(i - 1, j)]
            + grid.data[(i - 2, j)])
            / (twelve * grid.dx)
        4
