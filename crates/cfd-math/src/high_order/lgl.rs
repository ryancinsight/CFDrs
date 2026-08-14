//! Legendre-Gauss-Lobatto (LGL) nodes and weights — SSOT for `cfd-math` high-order methods.
//!
//! Gauss-Lobatto quadrature always includes the endpoints ±1.  An N-point rule is exact
//! for polynomials up to degree 2N − 3; the corresponding polynomial order is `N − 1`.
//!
//! Both the Discontinuous-Galerkin basis module and the Spectral-Element method module
//! require the same Newton-iteration LGL computation.  This module provides that single
//! implementation using `leto_ops::legendre_poly_and_deriv` as the Legendre-polynomial SSOT.
//!
//! Reference: Hesthaven & Warburton (2008) §A.3; Canuto et al. (2006) §A.1.

#![cfg_attr(
    test,
    expect(
        clippy::unwrap_used,
        reason = "ratchet CFDRS-UNWRAP-1: pre-existing debt"
    )
)]

use crate::error::Result;
use cfd_core::error::Error;
use leto_ops::{legendre_poly, legendre_poly_and_deriv};
use std::f64::consts::PI;

/// Compute `order + 1` Legendre-Gauss-Lobatto nodes on `[−1, 1]` for polynomial order `order`.
///
/// The two endpoints `±1` are fixed; the `order − 1` interior nodes are the roots of
/// `P_order′(x)`, found via Newton iteration with Chebyshev initial guesses.
pub fn lgl_nodes(order: usize) -> Result<Vec<f64>> {
    if order < 1 {
        return Err(Error::InvalidInput(format!(
            "LGL polynomial order must be at least 1, got {order}"
        )));
    }

    let n_pts = order + 1;
    let mut nodes = vec![0.0_f64; n_pts];
    nodes[0] = -1.0;
    nodes[order] = 1.0;

    // Chebyshev initial guesses for the interior roots of P_order′.
    for i in 1..order {
        nodes[i] = -((i as f64 * PI) / order as f64).cos();
    }

    // Newton iteration: find roots of P_order′ by solving P_order″(x) ΔΔx = P_order′(x).
    // P_order″(x) = (2x P_order′(x) − order(order+1) P_order(x)) / (1 − x²).
    const MAX_ITER: usize = 100;
    const TOL: f64 = 1e-15;

    for i in 1..order {
        let mut x = nodes[i];
        let mut delta = 1.0_f64;
        let mut iter = 0;

        while delta > TOL && iter < MAX_ITER {
            let (p, dp) = legendre_poly_and_deriv(order, x);
            let denom = 1.0 - x * x;
            if denom.abs() < 1e-15 {
                return Err(Error::Solver(format!(
                    "LGL Newton iteration hit endpoint at x = {x} for order = {order}"
                )));
            }
            let d2p = (2.0 * x * dp - (order * (order + 1)) as f64 * p) / denom;
            if d2p.abs() < f64::EPSILON {
                return Err(Error::Solver(format!(
                    "LGL Newton zero second derivative at x = {x} for order = {order}"
                )));
            }
            let dx = dp / d2p;
            x -= dx;
            delta = dx.abs();
            iter += 1;
        }

        if iter >= MAX_ITER {
            return Err(Error::Solver(format!(
                "LGL Newton failed to converge for interior node {i} of order {order}"
            )));
        }

        nodes[i] = x;
    }

    // Enforce symmetry to eliminate floating-point asymmetry from iteration.
    for i in 0..=(order / 2) {
        let val = nodes[i].abs();
        nodes[i] = -val;
        nodes[order - i] = val;
    }

    Ok(nodes)
}

/// Compute LGL quadrature weights for given nodes and polynomial order.
///
/// Uses the formula `w_i = 2 / (order(order+1) * P_order(x_i)²)`.
#[must_use]
pub fn lgl_weights(nodes: &[f64], order: usize) -> Vec<f64> {
    nodes
        .iter()
        .map(|&x| {
            let pn = legendre_poly(order, x);
            2.0 / (order as f64 * (order + 1) as f64) / (pn * pn)
        })
        .collect()
}

/// Compute `order + 1` LGL nodes *and* their quadrature weights in one call.
///
/// This is the primary entry point for both the DG-basis and spectral-element modules.
pub fn lgl_nodes_and_weights(order: usize) -> Result<(Vec<f64>, Vec<f64>)> {
    let nodes = lgl_nodes(order)?;
    let weights = lgl_weights(&nodes, order);
    Ok((nodes, weights))
}

#[cfg(test)]
mod tests {
    use super::*;
    use eunomia::assert_relative_eq;

    #[test]
    fn lgl_nodes_order2_gives_minus1_0_plus1() {
        let nodes = lgl_nodes(2).unwrap();
        assert_eq!(nodes.len(), 3);
        assert_relative_eq!(nodes[0], -1.0, epsilon = 1e-10);
        assert_relative_eq!(nodes[1], 0.0, epsilon = 1e-10);
        assert_relative_eq!(nodes[2], 1.0, epsilon = 1e-10);
    }

    #[test]
    fn lgl_nodes_order4_symmetric() {
        let nodes = lgl_nodes(4).unwrap();
        assert_eq!(nodes.len(), 5);
        assert_relative_eq!(nodes[0], -1.0, epsilon = 1e-10);
        assert_relative_eq!(nodes[4], 1.0, epsilon = 1e-10);
        assert_relative_eq!(nodes[1], -nodes[3], epsilon = 1e-10);
        assert_relative_eq!(nodes[2], 0.0, epsilon = 1e-10);
    }

    #[test]
    fn lgl_weights_order2_sum_to_two() {
        let nodes = lgl_nodes(2).unwrap();
        let weights = lgl_weights(&nodes, 2);
        assert_relative_eq!(weights[0], 1.0 / 3.0, epsilon = 1e-10);
        assert_relative_eq!(weights[1], 4.0 / 3.0, epsilon = 1e-10);
        assert_relative_eq!(weights[2], 1.0 / 3.0, epsilon = 1e-10);
        assert_relative_eq!(weights.iter().sum::<f64>(), 2.0, epsilon = 1e-10);
    }

    #[test]
    fn lgl_quadrature_exact_for_polynomials() {
        // An N-point LGL rule is exact for polynomials of degree ≤ 2N − 3.
        for order in 2usize..=5 {
            let (nodes, weights) = lgl_nodes_and_weights(order).unwrap();
            let n_pts = order + 1;
            // Maximum exact degree = 2*(n_pts) - 3 = 2*(order+1) - 3 = 2*order - 1.
            let max_exact_degree = 2 * order - 1;
            for d in 0..=max_exact_degree {
                let exact = if d % 2 == 0 {
                    2.0 / (d as f64 + 1.0)
                } else {
                    0.0
                };
                let integral: f64 = nodes
                    .iter()
                    .zip(weights.iter())
                    .map(|(&x, &w)| w * x.powi(d as i32))
                    .sum();
                assert!(
                    (integral - exact).abs() < 1e-10,
                    "LGL order={order} n_pts={n_pts} degree={d}: got {integral}, expected {exact}"
                );
            }
        }
    }
}
