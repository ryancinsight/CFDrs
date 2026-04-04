//! Gauss-Newton constraint solver for 2D parametric sketches.
//!
//! # Theorem — Gauss-Newton Convergence
//!
//! For a residual vector `r(x)` with Jacobian `J = dr/dx`, the
//! Gauss-Newton update is `dx = -(J^T J)^{-1} J^T r(x)`.
//! For a well-posed system where `J` has full column rank, the method
//! converges quadratically near the solution. For a fully constrained
//! sketch where `m = n`, `J` is square and this reduces to Newton's method.  ∎

use nalgebra::DMatrix;

use crate::application::sketch::residual::{constraint_residual_and_jacobian, JacobianRow};
use crate::domain::sketch::dof::{analyze_dofs, DofAnalysis};
use crate::domain::sketch::sketch::Sketch;

/// Gauss-Newton constraint solver.
pub struct ConstraintSolver {
    /// Maximum number of iterations.
    pub max_iterations: usize,
    /// Convergence tolerance on residual norm.
    pub tolerance: f64,
}

/// Result of a solve attempt.
#[derive(Clone, Debug)]
pub struct SolveResult {
    /// Whether the solver converged within tolerance.
    pub converged: bool,
    /// Number of iterations performed.
    pub iterations: usize,
    /// Final residual L2 norm.
    pub final_residual_norm: f64,
    /// DOF analysis after solving.
    pub dof_analysis: DofAnalysis,
}

impl Default for ConstraintSolver {
    fn default() -> Self {
        Self {
            max_iterations: 20,
            tolerance: 1e-10,
        }
    }
}

impl ConstraintSolver {
    /// Create a solver with default settings.
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    /// Solve the constraint system, modifying point positions in the sketch.
    pub fn solve(&self, sketch: &mut Sketch) -> SolveResult {
        let dof = analyze_dofs(sketch);
        if dof.total_dofs == 0 || dof.constraint_equations == 0 {
            return SolveResult {
                converged: true,
                iterations: 0,
                final_residual_norm: 0.0,
                dof_analysis: dof,
            };
        }

        let n_params = dof.total_dofs;
        let mut params = sketch.parameter_vector();

        let mut final_norm = f64::MAX;
        let mut iters = 0;

        for iteration in 0..self.max_iterations {
            iters = iteration + 1;

            let param_map = sketch.parameter_map();
            let constraints = sketch.constraints().to_vec();

            // Build residual and Jacobian.
            let mut residuals = Vec::new();
            let mut jac_rows: Vec<JacobianRow> = Vec::new();

            for (_, constraint) in &constraints {
                let (r, j) = constraint_residual_and_jacobian(constraint, sketch, &param_map);
                residuals.extend(r);
                jac_rows.extend(j);
            }

            let m = residuals.len();
            if m == 0 {
                final_norm = 0.0;
                break;
            }

            // Check convergence.
            final_norm = residuals.iter().map(|r| r * r).sum::<f64>().sqrt();
            if final_norm < self.tolerance {
                break;
            }

            // Assemble dense Jacobian matrix.
            let mut j_mat = DMatrix::zeros(m, n_params);
            for (row_idx, jrow) in jac_rows.iter().enumerate() {
                for &(col, val) in &jrow.entries {
                    if col < n_params {
                        j_mat[(row_idx, col)] = val;
                    }
                }
            }

            // Gauss-Newton: solve (J^T J) dx = -J^T r
            let jt = j_mat.transpose();
            let jtj = &jt * &j_mat;
            let r_vec = DMatrix::from_column_slice(m, 1, &residuals);
            let jtr = &jt * &r_vec;

            // Add small regularization for numerical stability.
            let reg = DMatrix::identity(n_params, n_params) * 1e-12;
            let lhs = jtj + reg;

            // Solve via LU decomposition.
            let lu = lhs.lu();
            let dx = lu.solve(&(-&jtr));
            let Some(dx) = dx else { break };

            // Update parameters with damping.
            let damping = 1.0;
            for i in 0..n_params {
                params[i] += damping * dx[(i, 0)];
            }

            sketch.apply_parameters(&params);
        }

        let final_dof = analyze_dofs(sketch);
        SolveResult {
            converged: final_norm < self.tolerance,
            iterations: iters,
            final_residual_norm: final_norm,
            dof_analysis: final_dof,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::domain::sketch::constraint::Constraint;
    use crate::domain::sketch::entity::{SketchEntity, SketchLine, SketchPoint};
    use crate::domain::sketch::sketch::{Sketch, SketchId};
    use crate::domain::sketch::work_plane::WorkPlane;
    use approx::assert_relative_eq;

    #[test]
    fn solver_converges_on_horizontal_constraint() {
        let mut sk = Sketch::new(SketchId(0), "test".into(), WorkPlane::xy());
        let p0 = sk.next_entity_id();
        let p1 = sk.next_entity_id();
        let ln = sk.next_entity_id();
        sk.add_entity(SketchEntity::Point(SketchPoint { id: p0, x: 0.0, y: 0.0, construction: false }));
        sk.add_entity(SketchEntity::Point(SketchPoint { id: p1, x: 3.0, y: 1.0, construction: false }));
        sk.add_entity(SketchEntity::Line(SketchLine { id: ln, start: p0, end: p1, construction: false }));
        sk.add_constraint(Constraint::Horizontal(ln));

        let solver = ConstraintSolver::new();
        let result = solver.solve(&mut sk);
        assert!(result.converged);

        let py0 = sk.point(p0).unwrap().y;
        let py1 = sk.point(p1).unwrap().y;
        assert_relative_eq!(py0, py1, epsilon = 1e-8);
    }

    #[test]
    fn solver_converges_on_distance_constraint() {
        let mut sk = Sketch::new(SketchId(0), "test".into(), WorkPlane::xy());
        let p0 = sk.next_entity_id();
        let p1 = sk.next_entity_id();
        sk.add_entity(SketchEntity::Point(SketchPoint { id: p0, x: 0.0, y: 0.0, construction: false }));
        sk.add_entity(SketchEntity::Point(SketchPoint { id: p1, x: 3.0, y: 0.0, construction: false }));
        sk.add_constraint(Constraint::Fixed(p0));
        sk.add_constraint(Constraint::Distance { a: p0, b: p1, value: 5.0, driving: true });

        let solver = ConstraintSolver::new();
        let result = solver.solve(&mut sk);
        assert!(result.converged);

        let px1 = sk.point(p1).unwrap().x;
        let py1 = sk.point(p1).unwrap().y;
        let dist = (px1 * px1 + py1 * py1).sqrt();
        assert_relative_eq!(dist, 5.0, epsilon = 1e-6);
    }

    #[test]
    fn empty_sketch_converges_immediately() {
        let mut sk = Sketch::new(SketchId(0), "empty".into(), WorkPlane::xy());
        let solver = ConstraintSolver::new();
        let result = solver.solve(&mut sk);
        assert!(result.converged);
        assert_eq!(result.iterations, 0);
    }
}
