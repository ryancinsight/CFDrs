use aequitas::systems::si::quantities::{
    Dimensionless, DynamicViscosity, Length, MassDensity, Pressure, SpecificHeatCapacity,
    TemperatureDifference, ThermalConductivity, ThermodynamicTemperature, Time, Velocity,
};

/// Plasma viscosity at 37°C [Pa·s]
/// Reference: Merrill et al. (1969)
pub const PLASMA_VISCOSITY_37C: DynamicViscosity = DynamicViscosity::from_base(0.00122);

/// Blood density in kg/m³.
///
/// ```compile_fail
/// use aequitas::systems::si::quantities::DynamicViscosity;
/// use cfd_core::physics::fluid::blood::constants::BLOOD_DENSITY;
/// let _: DynamicViscosity = BLOOD_DENSITY;
/// ```
/// Reference: Fung (1993)
pub const BLOOD_DENSITY: MassDensity = MassDensity::from_base(1060.0);

/// Zero-shear viscosity for normal blood (H_t = 45%) [Pa·s]
/// Reference: Cho & Kensey (1991)
pub const ZERO_SHEAR_VISCOSITY: DynamicViscosity = DynamicViscosity::from_base(0.056);

/// Infinite-shear viscosity for normal blood [Pa·s]
/// Reference: Cho & Kensey (1991)
pub const INFINITE_SHEAR_VISCOSITY: DynamicViscosity = DynamicViscosity::from_base(0.00345);

/// Yield stress for normal blood (H_t = 45%) \[Pa]
/// Reference: Merrill et al. (1969)
pub const YIELD_STRESS: Pressure = Pressure::from_base(0.0056);

/// Carreau-Yasuda relaxation time λ \[s]
/// Reference: Cho & Kensey (1991)
pub const CARREAU_LAMBDA: Time = Time::from_base(3.313);

/// Carreau-Yasuda power-law index n [-]
/// Reference: Cho & Kensey (1991)
pub const CARREAU_N: Dimensionless = Dimensionless::from_base(0.3568);

/// Carreau-Yasuda transition parameter a [-]
/// Reference: Cho & Kensey (1991)
pub const CARREAU_A: Dimensionless = Dimensionless::from_base(2.0);

/// Blood specific heat capacity [J/(kg·K)]
/// Reference: Fung (1993)
pub const BLOOD_SPECIFIC_HEAT: SpecificHeatCapacity = SpecificHeatCapacity::from_base(3770.0);

/// Blood thermal conductivity [W/(m·K)]
/// Reference: Fung (1993)
pub const BLOOD_THERMAL_CONDUCTIVITY: ThermalConductivity = ThermalConductivity::from_base(0.52);

/// Speed of sound in blood \[m/s]
/// Reference: Fung (1993)
pub const BLOOD_SPEED_OF_SOUND: Velocity = Velocity::from_base(1570.0);

/// Speed of sound in blood — medical ultrasound approximation \[m/s]
/// Commonly used value in therapeutic ultrasound literature (Duck 1990).
pub const BLOOD_SPEED_OF_SOUND_APPROX: Velocity = Velocity::from_base(1540.0);

/// Normal hematocrit (volume fraction of RBCs)
pub const NORMAL_HEMATOCRIT: Dimensionless = Dimensionless::from_base(0.45);

/// Critical vessel diameter for Fåhræus-Lindqvist effect \[m]
pub const FAHRAEUS_LINDQVIST_CRITICAL_DIAMETER: Length = Length::from_base(300e-6);

/// Normal human body temperature \[K] (37 C)
pub const BODY_TEMPERATURE: ThermodynamicTemperature = ThermodynamicTemperature::from_base(310.15);

/// Andrade blood viscosity temperature coefficient B \[K]
/// Used in the exponential temperature correction: exp(B * (1/T - 1/T_ref))
/// Reference: Fung (1993)
pub const ANDRADE_B_BLOOD: TemperatureDifference = TemperatureDifference::from_base(1500.0);
