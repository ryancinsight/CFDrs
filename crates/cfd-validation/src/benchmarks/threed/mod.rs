//! 3D benchmark submodules for CFD validation

pub mod bifurcation;
pub mod forced_turbulence;
pub mod nufft_coupling;
pub mod taylor_green;
pub mod serpentine;
pub mod venturi;

pub use bifurcation::BifurcationFlow3D;
pub use forced_turbulence::{
	ForcedTurbulenceBenchmark3D, ForcedTurbulenceBenchmarkConfig,
	ForcedTurbulenceBenchmarkHistory, ForcedTurbulenceBenchmarkReport,
};
pub use nufft_coupling::{
	NufftCouplingBenchmark3D, NufftCouplingBenchmarkConfig, NufftCouplingBenchmarkHistory,
	NufftCouplingBenchmarkReport,
};
pub use taylor_green::{
	TaylorGreenBenchmark3D, TaylorGreenBenchmarkConfig, TaylorGreenBenchmarkHistory,
	TaylorGreenBenchmarkReport,
};
pub use serpentine::SerpentineFlow3D;
pub use venturi::VenturiFlow3D;
