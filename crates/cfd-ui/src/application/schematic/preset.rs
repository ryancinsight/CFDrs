//! Preset loader — loads predefined millifluidic network topologies.

/// Available preset network topologies.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PresetType {
    SymmetricBifurcation,
    SymmetricTrifurcation,
    SerpentineChain,
    VenturiChain,
}

impl PresetType {
    /// Human-readable name for this preset.
    #[must_use]
    pub fn label(self) -> &'static str {
        match self {
            Self::SymmetricBifurcation => "Symmetric Bifurcation",
            Self::SymmetricTrifurcation => "Symmetric Trifurcation",
            Self::SerpentineChain => "Serpentine Chain",
            Self::VenturiChain => "Venturi Chain",
        }
    }

    /// All available presets.
    #[must_use]
    pub fn all() -> &'static [Self] {
        &[
            Self::SymmetricBifurcation,
            Self::SymmetricTrifurcation,
            Self::SerpentineChain,
            Self::VenturiChain,
        ]
    }
}
