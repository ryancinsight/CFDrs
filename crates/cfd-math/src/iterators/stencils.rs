//! Stencil operations for finite difference computations

use eunomia::{FloatElement, RealField};

fn from_f64<T: FloatElement>(value: f64) -> T {
    <T as FloatElement>::from_f64(value)
}

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
    pub fn first_derivative_coefficients<T: RealField + FloatElement>(&self) -> Vec<T> {
        match self {
            Self::Central3 => vec![from_f64(-0.5), T::ZERO, from_f64(0.5)],
            Self::Forward3 => vec![from_f64(-1.5), from_f64(2.0), from_f64(-0.5)],
            Self::Backward3 => vec![from_f64(0.5), from_f64(-2.0), from_f64(1.5)],
            Self::Central5 => vec![
                from_f64(1.0 / 12.0),
                from_f64(-2.0 / 3.0),
                T::ZERO,
                from_f64(2.0 / 3.0),
                from_f64(-1.0 / 12.0),
            ],
            Self::Upwind3 => vec![from_f64(-0.5), from_f64(-0.5), T::ONE],
        }
    }

    /// Get stencil coefficients for second derivative
    pub fn second_derivative_coefficients<T: RealField + FloatElement>(&self) -> Vec<T> {
        match self {
            Self::Central3 | Self::Forward3 | Self::Backward3 | Self::Upwind3 => {
                vec![T::ONE, from_f64(-2.0), T::ONE]
            }
            Self::Central5 => vec![
                from_f64(-1.0 / 12.0),
                from_f64(4.0 / 3.0),
                from_f64(-5.0 / 2.0),
                from_f64(4.0 / 3.0),
                from_f64(-1.0 / 12.0),
            ],
        }
    }

    /// Get the stencil size
    pub fn size(&self) -> usize {
        match self {
            Self::Central3 | Self::Forward3 | Self::Backward3 | Self::Upwind3 => 3,
            Self::Central5 => 5,
        }
    }
}

/// Iterator for applying stencils to data
pub struct StencilIterator<I, T> {
    iter: I,
    pattern: StencilPattern,
    buffer: Vec<T>,
}

impl<I, T> StencilIterator<I, T>
where
    I: Iterator<Item = T>,
    T: RealField + FloatElement + Copy,
{
    /// Create a new stencil iterator
    pub fn new(iter: I, pattern: StencilPattern) -> Self {
        Self {
            iter,
            pattern,
            buffer: Vec::with_capacity(pattern.size()),
        }
    }

    /// Apply stencil for first derivative
    pub fn first_derivative(&mut self) -> Option<T> {
        self.apply_stencil(self.pattern.first_derivative_coefficients())
    }

    /// Apply stencil for second derivative
    pub fn second_derivative(&mut self) -> Option<T> {
        self.apply_stencil(self.pattern.second_derivative_coefficients())
    }

    fn apply_stencil(&mut self, coefficients: Vec<T>) -> Option<T> {
        // Fill buffer if needed
        while self.buffer.len() < self.pattern.size() {
            let value = self.iter.next()?;
            self.buffer.push(value);
        }

        // Apply coefficients
        let result = self
            .buffer
            .iter()
            .zip(coefficients.iter())
            .map(|(val, coeff)| *val * *coeff)
            .fold(T::ZERO, |acc, x| acc + x);

        // Slide buffer
        if let Some(next) = self.iter.next() {
            self.buffer.remove(0);
            self.buffer.push(next);
        }

        Some(result)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use eunomia::assert_relative_eq;

    #[test]
    fn first_derivative_coefficients_match_known_stencils() {
        assert_eq!(
            StencilPattern::Forward3.first_derivative_coefficients::<f64>(),
            vec![-1.5, 2.0, -0.5]
        );
        assert_eq!(
            StencilPattern::Backward3.first_derivative_coefficients::<f64>(),
            vec![0.5, -2.0, 1.5]
        );

        let central5 = StencilPattern::Central5.first_derivative_coefficients::<f64>();
        assert_relative_eq!(central5[0], 1.0 / 12.0, epsilon = 1e-15);
        assert_relative_eq!(central5[1], -2.0 / 3.0, epsilon = 1e-15);
        assert_relative_eq!(central5[3], 2.0 / 3.0, epsilon = 1e-15);
        assert_relative_eq!(central5[4], -1.0 / 12.0, epsilon = 1e-15);
    }

    #[test]
    fn second_derivative_coefficients_are_input_sensitive() {
        let expected = vec![1.0, -2.0, 1.0];
        assert_eq!(
            StencilPattern::Central3.second_derivative_coefficients::<f64>(),
            expected
        );
        assert_eq!(
            StencilPattern::Forward3.second_derivative_coefficients::<f64>(),
            expected
        );
        assert_eq!(
            StencilPattern::Backward3.second_derivative_coefficients::<f64>(),
            expected
        );
        assert_eq!(
            StencilPattern::Upwind3.second_derivative_coefficients::<f64>(),
            expected
        );
    }

    #[test]
    fn stencil_iterator_applies_coefficients() {
        let values = vec![1.0, 2.0, 4.0];
        let mut iter = StencilIterator::new(values.into_iter(), StencilPattern::Forward3);

        assert_relative_eq!(iter.first_derivative().unwrap(), 0.5, epsilon = 1e-15);
    }
}
