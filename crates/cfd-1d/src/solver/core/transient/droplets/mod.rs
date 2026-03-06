//! Transient droplet and multiphase flow tracking
//!
//! ## Theorem: Discrete 1D Multiphase Droplet Tracking
//!
//! **Theorem**: A dispersed phase (droplets) in a continuous 1D carrier fluid
//! transports kinematically according to the superficial fluid velocity $u(t)$,
//! modified by a slip factor $\beta$ due to wall lubrication layer effects:
//!
//! $$ \frac{dx_d}{dt} = \beta u(t) = \beta \frac{Q(t)}{A(x)} $$
//!
//! **Conservation of Droplet Volume**:
//! The volume $V_d$ of an incompressible droplet is an invariant scalar across
//! topology changes. At a bifurcation junction where the superficial flow $Q$
//! splits into $Q_1$ and $Q_2$, a droplet undergoes deterministic fragmentation
//! proportional to the time-integrated flow split if the capillary number exceeds
//! the critical threshold for breakup:
//!
//! $$ V_{d,1} = V_d \left( \frac{Q_1}{Q_1 + Q_2} \right), \quad V_{d,2} = V_d \left( \frac{Q_2}{Q_1 + Q_2} \right) $$

/// Main droplet simulation engine
pub mod simulator;
/// Droplet and split policy types
pub mod types;

pub use simulator::TransientDropletSimulator;
pub use types::{
    ChannelOccupancy, DropletBoundary, DropletInjection, DropletPosition, DropletSnapshot,
    DropletSplitPolicy, DropletState, DropletTrackingState, SplitMode,
};
