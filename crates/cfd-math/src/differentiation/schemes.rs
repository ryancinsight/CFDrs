//! Finite difference schemes for numerical differentiation.

/// Finite difference schemes
#[derive(Debug, Clone, Copy, PartialEq)]
#[derive(Default)]
pub enum FiniteDifferenceScheme {
    /// Forward difference: f'(x) ≈ (f(x+h) - f(x)) / h
    Forward,
    /// Backward difference: f'(x) ≈ (f(x) - f(x-h)) / h
    Backward,
    /// Central difference: f'(x) ≈ (f(x+h) - f(x-h)) / (2h)
    #[default]
    Central,
    /// Second-order forward: f'(x) ≈ (-3f(x) + 4f(x+h) - f(x+2h)) / (2h)
    ForwardSecondOrder,
    /// Second-order backward: f'(x) ≈ (f(x-2h) - 4f(x-h) + 3f(x)) / (2h)
    BackwardSecondOrder,
}

