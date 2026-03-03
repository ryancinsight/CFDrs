//! `PyO3` Python bindings for CFD-rs solvers
//!
//! This crate exposes the Rust CFD solvers to Python for easy validation,
//! comparison with other CFD codes, and interactive analysis.
//!
//! # Building
//!
//! ```bash
//! maturin develop
//! ```
//!
//! # Python Usage
//!
//! ```python
//! import cfd_python
//!
//! # Create bifurcation solver
//! bifurc = cfd_python.BifurcationSolver(d_parent=100e-6, d_daughter1=80e-6, d_daughter2=80e-6)
//!
//! # Create blood model
//! blood = cfd_python.CassonBlood()
//!
//! # Solve
//! result = bifurc.solve(flow_rate=3e-8, pressure=40.0, blood=blood)
//! print(f"Pressure drop 1: {result.dp_1} Pa")
//! print(f"Wall shear rate 1: {result.gamma_1} s^-1")
//! ```

#![warn(clippy::all)]
#![warn(clippy::pedantic)]
// PyO3 bindings require specific patterns that trigger pedantic lints.
#![allow(clippy::cast_precision_loss)]
#![allow(clippy::cast_possible_truncation)]
#![allow(clippy::cast_sign_loss)]
#![allow(clippy::similar_names)]
#![allow(clippy::many_single_char_names)]
#![allow(clippy::too_many_arguments)]
#![allow(clippy::too_many_lines)]
#![allow(clippy::module_name_repetitions)]
#![allow(clippy::missing_errors_doc)]
#![allow(clippy::missing_panics_doc)]
#![allow(clippy::must_use_candidate)]
#![allow(clippy::needless_pass_by_value)]
#![allow(clippy::redundant_closure_for_method_calls)]
#![allow(clippy::doc_markdown)]
#![allow(clippy::unreadable_literal)]
#![allow(clippy::used_underscore_binding)]
// PyO3 getters/setters don't follow Rust naming conventions.
#![allow(clippy::unused_self)]
#![allow(clippy::manual_let_else)]
#![allow(clippy::unnecessary_wraps)]
#![allow(clippy::match_same_arms)]
#![allow(clippy::useless_conversion)]
#![allow(clippy::inline_always)]
#![allow(clippy::needless_range_loop)]
#![allow(clippy::implicit_hasher)]
#![allow(clippy::new_without_default)]
#![allow(clippy::trivially_copy_pass_by_ref)]
#![allow(clippy::ptr_arg)]
#![allow(clippy::format_push_string)]
#![allow(clippy::field_reassign_with_default)]
#![allow(clippy::empty_line_after_doc_comments)]
#![allow(clippy::len_without_is_empty)]
#![allow(clippy::new_ret_no_self)]
#![allow(clippy::items_after_statements)]
#![allow(clippy::float_cmp)]
#![allow(clippy::return_self_not_must_use)]
#![allow(clippy::cast_possible_wrap)]
#![allow(clippy::doc_overindented_list_items)]
#![allow(clippy::same_item_push)]

use pyo3::prelude::*;

mod bifurcation;
mod blood;
mod poiseuille_2d;
mod result_types;
mod solver_2d;
mod solver_3d;
mod womersley;

pub use bifurcation::{PyBifurcationSolver, PyTrifurcationResult, PyTrifurcationSolver};
pub use blood::*;
pub use womersley::*;
pub use poiseuille_2d::{PyPoiseuilleConfig, PyPoiseuilleResult, PyPoiseuilleSolver};
pub use result_types::PyBifurcationResult;
pub use solver_2d::*;
pub use solver_3d::*;

/// `PyO3` module for CFD-rs Python bindings
#[pymodule]
fn cfd_python(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add("__version__", "0.1.0")?;

    // 1D solvers
    m.add_class::<PyBifurcationSolver>()?;
    m.add_class::<PyBifurcationResult>()?;
    m.add_class::<PyTrifurcationSolver>()?;
    m.add_class::<PyTrifurcationResult>()?;

    // 2D solvers
    m.add_class::<PyPoiseuilleSolver>()?;
    m.add_class::<PyPoiseuilleResult>()?;
    m.add_class::<PyPoiseuilleConfig>()?;

    // Blood models
    m.add_class::<PyCassonBlood>()?;
    m.add_class::<PyCarreauYasudaBlood>()?;
    m.add_class::<PyCrossBlood>()?;
    m.add_class::<PyFahraeuasLindqvist>()?;

    // Womersley pulsatile flow
    m.add_class::<PyWomersleyNumber>()?;
    m.add_class::<PyWomersleyProfile>()?;
    m.add_class::<PyWomersleyFlow>()?;

    // 2D solvers (extended)
    m.add_class::<PyPoiseuille2DSolver>()?;
    m.add_class::<PyPoiseuille2DResult>()?;
    m.add_class::<PyVenturiSolver2D>()?;
    m.add_class::<PyVenturiResult2D>()?;
    m.add_class::<PyTrifurcationSolver2D>()?;
    m.add_class::<PyTrifurcationResult2D>()?;

    // 3D solvers
    m.add_class::<PyBifurcation3DSolver>()?;
    m.add_class::<PyBifurcation3DResult>()?;
    m.add_class::<PyTrifurcation3DSolver>()?;
    m.add_class::<PyTrifurcation3DResult>()?;
    m.add_class::<PyPoiseuille3DSolver>()?;
    m.add_class::<PyPoiseuille3DResult>()?;
    m.add_class::<PyVenturi3DSolver>()?;
    m.add_class::<PyVenturi3DResult>()?;
    m.add_class::<PySerpentine3DSolver>()?;
    m.add_class::<PySerpentine3DResult>()?;

    // 2D Cavity solver (Ghia benchmark)
    m.add_class::<PyCavitySolver2D>()?;
    m.add_class::<PyCavityResult2D>()?;

    // 2D Bifurcation solver
    m.add_class::<PyBifurcationSolver2D>()?;
    m.add_class::<PyBifurcationResult2D>()?;

    // 1D Serpentine resistance solver
    m.add_class::<PySerpentineSolver1D>()?;
    m.add_class::<PySerpentineResult1D>()?;

    // 1D Venturi resistance solver
    m.add_class::<PyVenturiSolver1D>()?;
    m.add_class::<PyVenturiResult1D>()?;

    // Add submodule for validation utilities
    let validation = PyModule::new_bound(m.py(), "validation")?;
    validation.add("__doc__", "Validation utilities for CFD comparisons")?;
    m.add_submodule(&validation)?;

    Ok(())
}
