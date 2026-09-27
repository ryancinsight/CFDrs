/// DES model variants
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum DESVariant {
    /// Original DES97
    DES97,
    /// Delayed DES (DDES)
    DDES,
    /// Improved DDES (IDDES)
    IDDES,
}

/// DES configuration
#[derive(Debug, Clone)]
pub struct DESConfig {
    /// DES variant to use
    pub variant: DESVariant,
    /// DES constant (C_DES)
    pub des_constant: f64,
    /// Maximum SGS viscosity ratio
    pub max_sgs_ratio: f64,
    /// Molecular viscosity (constant)
    pub rans_viscosity: f64,
    /// Enable GPU acceleration
    pub use_gpu: bool,
}

impl Default for DESConfig {
    fn default() -> Self {
        Self {
            variant: DESVariant::DDES,
            des_constant: 0.65,
            max_sgs_ratio: 0.5,   // Prevent excessive SGS viscosity
            rans_viscosity: 1e-5, // Default molecular viscosity
            use_gpu: false,
        }
    }
}
