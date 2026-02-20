pub mod types;
pub mod simulator;

pub use simulator::TransientDropletSimulator;
pub use types::{
    ChannelOccupancy, DropletBoundary, DropletInjection, DropletPosition, DropletSnapshot,
    DropletSplitPolicy, DropletState, DropletTrackingState, SplitMode,
};
