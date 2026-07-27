//! Transient composition and concentration tracking
//!
//! ## Theorem: Multispecies 1D Advective Mass Transport
//!
//! **Theorem**: In a 1D channel network where longitudinal diffusion is negligible
//! compared to advection ($\text{Pe} \gg 1$), the transport of non-reacting species
//! concentration $C(x,t)$ obeys the pure advection equation:
//!
//! $$ \frac{\partial C}{\partial t} + u(t) \frac{\partial C}{\partial x} = 0 $$
//!
//! **Discrete State Theorem**: For a well-mixed uniform segment control volume $V_e$,
//! the rate of change of species mass $M_e = C_e V_e$ strictly satisfies the
//! discrete mass balance:
//!
//! $$ \frac{d M_e}{dt} = \sum_{in} Q_{in} C_{in} - \sum_{out} Q_{out} C_e $$
//!
//! Exact volumetric tracking is achieved by discretizing time such that the Courant
//! number $\text{CFL} = \frac{u \Delta t}{L} \le 1$, ensuring no mass overshoots
//! the volume boundaries during a single timestep.

/// Events for flow and boundary changes.
pub mod events;
/// Main composition simulation engine.
pub mod simulator;
/// Mixture state definitions.
pub mod state;

pub use events::{
    EdgeFlowEvent, InletCompositionEvent, InletHematocritEvent, PressureBoundaryEvent,
};
pub use simulator::{
    BloodEdgeTransportConfig, SimulationTimeConfig, TransientCompositionSimulator,
};
pub use state::{CompositionState, MixtureComposition, BLOOD_PLASMA_FLUID_ID, BLOOD_RBC_FLUID_ID};
