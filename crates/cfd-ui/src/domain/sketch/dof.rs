//! Degree-of-freedom analysis for sketch constraint systems.

use super::sketch::Sketch;

/// Result of DOF analysis on a sketch.
#[derive(Clone, Debug)]
pub struct DofAnalysis {
    /// Total free parameters (2 per point + 1 per circle radius).
    pub total_dofs: usize,
    /// Total constraint equations.
    pub constraint_equations: usize,
    /// Overall constraint status.
    pub status: DofStatus,
    /// Remaining free DOFs (total_dofs - constraint_equations).
    pub free_dof_count: i32,
}

/// Constraint status of a sketch.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum DofStatus {
    /// More DOFs than constraints — geometry can still move.
    UnderConstrained,
    /// DOFs equal constraints — geometry is fully locked.
    FullyConstrained,
    /// More constraints than DOFs — conflicting constraints.
    OverConstrained,
}

/// Analyze the degrees of freedom of a sketch.
#[must_use]
pub fn analyze_dofs(sketch: &Sketch) -> DofAnalysis {
    let total_dofs = sketch.parameter_vector().len();
    let constraint_equations: usize = sketch
        .constraints()
        .iter()
        .map(|(_, c)| c.residual_count())
        .sum();

    let free = total_dofs as i32 - constraint_equations as i32;
    let status = if free > 0 {
        DofStatus::UnderConstrained
    } else if free == 0 {
        DofStatus::FullyConstrained
    } else {
        DofStatus::OverConstrained
    };

    DofAnalysis {
        total_dofs,
        constraint_equations,
        status,
        free_dof_count: free,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::domain::sketch::constraint::Constraint;
    use crate::domain::sketch::entity::{EntityId, SketchEntity, SketchPoint};
    use crate::domain::sketch::sketch::{Sketch, SketchId};
    use crate::domain::sketch::work_plane::WorkPlane;

    #[test]
    fn empty_sketch_has_zero_dofs() {
        let sk = Sketch::new(SketchId(0), "empty".into(), WorkPlane::xy());
        let analysis = analyze_dofs(&sk);
        assert_eq!(analysis.total_dofs, 0);
        assert_eq!(analysis.free_dof_count, 0);
        assert_eq!(analysis.status, DofStatus::FullyConstrained);
    }

    #[test]
    fn single_point_is_under_constrained() {
        let mut sk = Sketch::new(SketchId(0), "test".into(), WorkPlane::xy());
        let id = sk.next_entity_id();
        sk.add_entity(SketchEntity::Point(SketchPoint {
            id, x: 0.0, y: 0.0, construction: false,
        }));
        let analysis = analyze_dofs(&sk);
        assert_eq!(analysis.total_dofs, 2);
        assert_eq!(analysis.free_dof_count, 2);
        assert_eq!(analysis.status, DofStatus::UnderConstrained);
    }

    #[test]
    fn fixed_point_is_fully_constrained() {
        let mut sk = Sketch::new(SketchId(0), "test".into(), WorkPlane::xy());
        let id = sk.next_entity_id();
        sk.add_entity(SketchEntity::Point(SketchPoint {
            id, x: 0.0, y: 0.0, construction: false,
        }));
        sk.add_constraint(Constraint::Fixed(id));
        let analysis = analyze_dofs(&sk);
        assert_eq!(analysis.status, DofStatus::FullyConstrained);
    }
}
