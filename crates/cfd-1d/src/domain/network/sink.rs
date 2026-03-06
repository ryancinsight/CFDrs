//! [`GraphSink`] adapter: converts a [`NetworkBlueprint`] into a solver-ready
//! [`Network<T, F>`] using the Hagen-Poiseuille / Shah-London resistance model.
//!
//! This is the **canonical external entry-point** for consuming a blueprint from
//! `cfd-schematics` in a 1D Kirchhoff simulation.  The adapter pattern keeps all
//! blueprint *construction* logic inside `cfd-schematics` while `cfd-1d` only
//! needs to know how to convert a finished blueprint into its own graph type.
//!
//! # Usage
//! ```rust,ignore
//! use cfd_1d::NetworkBuilderSink;
//! use cfd_schematics::application::use_cases::generate_network::NetworkGenerationService;
//! use cfd_schematics::interface::presets::venturi_rect;
//!
//! let bp = venturi_rect("v", 2e-3, 0.5e-3, 0.5e-3, 2e-3);
//! let sink = NetworkBuilderSink::new(my_fluid);
//! let network = NetworkGenerationService::new(sink).generate(&bp)?;
//! ```

use std::marker::PhantomData;

use cfd_core::error::Result;
use cfd_core::physics::fluid::FluidTrait;
use cfd_schematics::application::ports::GraphSink;
use cfd_schematics::domain::model::NetworkBlueprint;
use nalgebra::RealField;
use num_traits::FromPrimitive;

use super::builder::network_from_blueprint;
use crate::domain::network::wrapper::Network;

/// A [`GraphSink`] that converts a validated [`NetworkBlueprint`] into a
/// solver-ready [`Network<T, F>`] using Hagen-Poiseuille / Shah-London
/// resistance models.
///
/// Construct with `NetworkBuilderSink::new(fluid)` and pass to
/// [`NetworkGenerationService`](cfd_schematics::application::use_cases::generate_network::NetworkGenerationService).
pub struct NetworkBuilderSink<T, F> {
    fluid: F,
    _t: PhantomData<T>,
}

impl<T, F> NetworkBuilderSink<T, F> {
    /// Create a new sink with the given fluid model.
    #[must_use]
    pub fn new(fluid: F) -> Self {
        Self {
            fluid,
            _t: PhantomData,
        }
    }
}

impl<T, F> GraphSink for NetworkBuilderSink<T, F>
where
    T: RealField + Copy + FromPrimitive,
    F: FluidTrait<T> + Clone,
{
    type Output = Network<T, F>;

    fn build(&self, blueprint: &NetworkBlueprint) -> Result<Network<T, F>> {
        network_from_blueprint(blueprint, self.fluid.clone())
    }
}
