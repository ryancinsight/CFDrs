//! Stencil operations for finite difference computations

use nalgebra::RealField;
/// Stencil patterns for finite difference schemes
#[derive(Debug, Clone, Copy)]
pub enum StencilPattern {
    /// Central difference (3-point)
    Central3,
    /// Forward difference (3-point)
    Forward3,
    /// Backward difference (3-point)
    Backward3,
    /// Central difference (5-point)
    Central5,
    /// Upwind (3-point)
    Upwind3,
}
impl StencilPattern {
    /// Get stencil coefficients for first derivative
    pub fn first_derivative_coefficients<T: RealField>(&self) -> Vec<T> {
        match self {
            Self::Central3 => vec![
                T::from_f64(-0.5).unwrap_or_else(T::zero),
                T::zero(),
                T::from_f64(0.5).unwrap_or_else(T::zero),
            ],
            Self::Forward3 => vec![
                T::from_f64(-1.5).unwrap_or_else(T::zero),
                T::from_f64(2.0).unwrap_or_else(T::zero),
            Self::Backward3 => vec![
                T::from_f64(-2.0).unwrap_or_else(T::zero),
                T::from_f64(1.5).unwrap_or_else(T::zero),
            Self::Central5 => vec![
                T::from_f64(1.0 / 12.0).unwrap_or_else(T::zero),
                T::from_f64(-2.0 / 3.0).unwrap_or_else(T::zero),
                T::from_f64(2.0 / 3.0).unwrap_or_else(T::zero),
                T::from_f64(-1.0 / 12.0).unwrap_or_else(T::zero),
            Self::Upwind3 => vec![
                T::one(),
        }
    }
    /// Get stencil coefficients for second derivative
    pub fn second_derivative_coefficients<T: RealField>(&self) -> Vec<T> {
                T::from_f64(4.0 / 3.0).unwrap_or_else(T::zero),
                T::from_f64(-5.0 / 2.0).unwrap_or_else(T::zero),
            _ => vec![T::zero(); 3], // Not implemented for other patterns
    /// Get the stencil size
    }

    pub fn size(&self) -> usize {
            Self::Central3 | Self::Forward3 | Self::Backward3 | Self::Upwind3 => 3,
            Self::Central5 => 5,
/// Iterator for applying stencils to data
    }

}

pub struct StencilIterator<I, T> {
    iter: I,
    pattern: StencilPattern,
    buffer: Vec<T>,
impl<I, T> StencilIterator<I, T>
where
    I: Iterator<Item = T>,
    T: RealField + Copy,
{
    /// Create a new stencil iterator
    pub fn new(iter: I, pattern: StencilPattern) -> Self {
        Self {
            iter,
            pattern,
            buffer: Vec::with_capacity(pattern.size()),
    /// Apply stencil for first derivative
    }

    pub fn first_derivative(&mut self) -> Option<T> {
        self.apply_stencil(self.pattern.first_derivative_coefficients())
    /// Apply stencil for second derivative
    }

    pub fn second_derivative(&mut self) -> Option<T> {
        self.apply_stencil(self.pattern.second_derivative_coefficients())
    }

    fn apply_stencil(&mut self, coefficients: Vec<T>) -> Option<T> {
        // Fill buffer if needed
        while self.buffer.len() < self.pattern.size() {
            if let Some(val) = self.iter.next() {
                self.buffer.push(val);
            } else {
                return None;
            }
        // Apply coefficients
        let result = self
            .buffer
            .iter()
            .zip(coefficients.iter())
            .map(|(val, coeff)| *val * *coeff)
            .fold(T::zero(), |acc, x| acc + x);
        // Slide buffer
        if let Some(next) = self.iter.next() {
            self.buffer.remove(0);
            self.buffer.push(next);
        Some(result)


}
}
}
}
}
}
}
}
