/// Plasma viscosity at 37°C [Pa·s]
/// Reference: Merrill et al. (1969)
pub const PLASMA_VISCOSITY_37C: f64 = 0.00122;

/// Blood density [kg/m³]
/// Reference: Fung (1993)
pub const BLOOD_DENSITY: f64 = 1060.0;

/// Zero-shear viscosity for normal blood (H_t = 45%) [Pa·s]
/// Reference: Cho & Kensey (1991)
pub const ZERO_SHEAR_VISCOSITY: f64 = 0.056;

/// Infinite-shear viscosity for normal blood [Pa·s]
/// Reference: Cho & Kensey (1991)
pub const INFINITE_SHEAR_VISCOSITY: f64 = 0.00345;

/// Yield stress for normal blood (H_t = 45%) [Pa]
/// Reference: Merrill et al. (1969)
pub const YIELD_STRESS: f64 = 0.0056;

/// Casson viscosity parameter √μ_∞ for normal blood [√(Pa·s)]
/// Reference: Merrill et al. (1969)
pub const CASSON_VISCOSITY_SQRT: f64 = 0.0588; // √0.00345 ≈ 0.0587

/// Carreau-Yasuda relaxation time λ [s]
/// Reference: Cho & Kensey (1991)
pub const CARREAU_LAMBDA: f64 = 3.313;

/// Carreau-Yasuda power-law index n [-]
/// Reference: Cho & Kensey (1991)
pub const CARREAU_N: f64 = 0.3568;

/// Carreau-Yasuda transition parameter a [-]
/// Reference: Cho & Kensey (1991)
pub const CARREAU_A: f64 = 2.0;

/// Blood specific heat capacity [J/(kg·K)]
/// Reference: Fung (1993)
pub const BLOOD_SPECIFIC_HEAT: f64 = 3770.0;

/// Blood thermal conductivity [W/(m·K)]
/// Reference: Fung (1993)
pub const BLOOD_THERMAL_CONDUCTIVITY: f64 = 0.52;

/// Speed of sound in blood [m/s]
/// Reference: Fung (1993)
pub const BLOOD_SPEED_OF_SOUND: f64 = 1570.0;

/// Normal hematocrit (volume fraction of RBCs)
pub const NORMAL_HEMATOCRIT: f64 = 0.45;

/// Critical vessel diameter for Fåhræus-Lindqvist effect [m]
pub const FAHRAEUS_LINDQVIST_CRITICAL_DIAMETER: f64 = 300e-6;
