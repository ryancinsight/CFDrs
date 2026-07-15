//! 2D CFD solver `PyO3` wrappers.
//!
//! This module exposes 2D finite volume solvers and 1D resistance models
//! for validation against Python CFD packages.
//!
//! ## Sub-modules
//!
//! | Module | Solvers |
//! |--------|---------|
//! | [`poiseuille`] | `Poiseuille2DSolver`, `Poiseuille2DResult` |
//! | [`venturi`] | `VenturiSolver2D/1D`, `VenturiResult2D/1D` |
//! | [`cavity`] | `CavitySolver2D`, `CavityResult2D` |
//! | [`bifurcation`] | `BifurcationSolver2D`, `BifurcationResult2D`, `TrifurcationSolver2D`, `TrifurcationResult2D` |
//! | [`serpentine`] | `SerpentineSolver1D`, `SerpentineResult1D` |

use leto::Array2;
use numpy::PyArray2;
use pyo3::exceptions::PyRuntimeError;
use pyo3::prelude::*;

mod bifurcation;
mod cavity;
mod poiseuille;
mod serpentine;
mod venturi;

pub use bifurcation::*;
pub use cavity::*;
pub use poiseuille::*;
pub use serpentine::*;
pub use venturi::*;

pub(super) fn pyarray2_from_leto(
    py: Python<'_>,
    array: Array2<f64>,
) -> PyResult<Bound<'_, PyArray2<f64>>> {
    let [rows, cols] = array.shape();
    if cols == 0 && rows != 0 {
        return Err(PyRuntimeError::new_err(
            "cannot expose non-empty zero-column Leto array as NumPy vec2",
        ));
    }

    let values = array.into_vec();
    let rows: Vec<Vec<f64>> = if cols == 0 {
        Vec::new()
    } else {
        values.chunks(cols).map(<[f64]>::to_vec).collect()
    };

    PyArray2::from_vec2_bound(py, &rows).map_err(|error| {
        PyRuntimeError::new_err(format!("NumPy array construction failed: {error}"))
    })
}
