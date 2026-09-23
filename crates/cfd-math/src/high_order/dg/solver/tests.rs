use super::super::DGOperator;
use crate::error::Result;
use leto::Array2;

use super::super::{matrix_from_element, matrix_neg, vector_from_vec};
use super::*;
use eunomia::assert_relative_eq;

use crate::high_order::{DGOperatorParams, FluxType, LimiterType};
type Jacobian = dyn Fn(f64, &Array2<f64>) -> Result<Array2<f64>>;

#[test]
fn test_forward_euler() {
    let f = |_t: f64, u: &Array2<f64>| Ok(matrix_neg(u));

    let dt = 0.01;
    let t_final = 1.0;
    let n_steps = (t_final / dt) as usize;

    let mut u = matrix_from_element(1, 1, 1.0);

    let integrator = ForwardEuler;

    let jac: Option<&Jacobian> = None;
    for _ in 0..n_steps {
        let (u_new, _) = integrator
            .step(0.0, dt, &u, &f, jac)
            .expect("expected value");
        u = u_new;
    }

    let exact = (-t_final).exp();
    assert_relative_eq!(u[[0, 0]], exact, epsilon = 1e-2);
}

#[test]
fn test_rk4() {
    let f = |_t: f64, u: &Array2<f64>| Ok(matrix_neg(u));

    let dt = 0.1;
    let t_final = 1.0;
    let n_steps = (t_final / dt) as usize;

    let mut u = matrix_from_element(1, 1, 1.0);

    let integrator = RK4;

    let jac: Option<&Jacobian> = None;
    for _ in 0..n_steps {
        let (u_new, _) = integrator
            .step(0.0, dt, &u, &f, jac)
            .expect("expected value");
        u = u_new;
    }

    let exact = (-t_final).exp();
    assert_relative_eq!(u[[0, 0]], exact, epsilon = 1e-6);
}

#[test]
fn test_ssprk3() {
    let f = |_t: f64, u: &Array2<f64>| Ok(matrix_neg(u));

    let dt = 0.1;
    let t_final = 1.0;
    let n_steps = (t_final / dt) as usize;

    let mut u = matrix_from_element(1, 1, 1.0);

    let integrator = SSPRK3;

    let jac: Option<&Jacobian> = None;
    for _ in 0..n_steps {
        let (u_new, _) = integrator
            .step(0.0, dt, &u, &f, jac)
            .expect("expected value");
        u = u_new;
    }

    let exact = (-t_final).exp();
    assert_relative_eq!(u[[0, 0]], exact, epsilon = 1e-4);
}

#[test]
fn test_implicit_euler() {
    let f = |_t: f64, u: &Array2<f64>| Ok(matrix_neg(u));
    let jac = |_t: f64, _u: &Array2<f64>| Ok(matrix_from_element(1, 1, -1.0));

    let dt = 0.1;
    let t_final = 1.0;
    let n_steps = (t_final / dt) as usize;

    let mut u = matrix_from_element(1, 1, 1.0);

    let integrator = ImplicitEuler::default();

    for _ in 0..n_steps {
        let (u_new, _) = integrator
            .step(0.0, dt, &u, &f, Some(&jac))
            .expect("expected value");
        u = u_new;
    }

    let exact = (-t_final).exp();
    assert_relative_eq!(u[[0, 0]], exact, epsilon = 2e-2);
}

#[test]
fn test_dg_solver() -> Result<()> {
    let order = 2;
    let num_components = 1;
    let params = DGOperatorParams::new()
        .with_volume_flux(FluxType::Central)
        .with_surface_flux(FluxType::LaxFriedrichs)
        .with_limiter(LimiterType::None);

    let dg_op = DGOperator::new(order, num_components, Some(params))?;

    let integrator = TimeIntegratorFactory::create(TimeIntegration::SSPRK3);

    let t_final = 1.0;
    let solver_params = TimeIntegrationParams::new(TimeIntegration::SSPRK3)
        .with_t_final(t_final)
        .with_dt(0.01)
        .with_verbose(false);

    let mut solver = DGSolver::new(dg_op, integrator, solver_params);

    let u0 = |x: f64| vector_from_vec(vec![1.0 + x + x * x]);
    solver.initialize(u0)?;

    let f = |_t: f64, u: &Array2<f64>| Ok(matrix_neg(u));
    solver.solve(f, None::<fn(f64, &Array2<f64>) -> Result<Array2<f64>>>)?;

    let x = 0.5;
    let u_num = solver.evaluate(x)[0];
    let u_exact = (1.0 + x + x * x) * (-t_final).exp();

    assert_relative_eq!(u_num, u_exact, epsilon = 1e-3);

    Ok(())
}
