//! Acoustic radiation force on blood cells (Gor'kov 1962).
//!
//! ## Theorem — Gor'kov Acoustic Radiation Force
//!
//! A small sphere (cell) of radius a in an acoustic field experiences
//! a time-averaged radiation force:
//!
//! ```text
//! F_rad = −(4π/3) a³ [β_m/(2ρ_m c_m²) · ∇⟨p²⟩ − (3ρ_m/2) · f₁ · ∇⟨v²⟩]
//! ```
//!
//! For the primary radiation force in a standing wave (frequency f, wavelength λ):
//!
//! ```text
//! F_rad = 4π a³ Φ k E_ac sin(2kx)
//! ```
//!
//! where:
//! - k = 2π/λ = 2πf/c (wavenumber)
//! - E_ac = p₀²/(4ρc²) (acoustic energy density) \[J/m³\]
//! - Φ = f₁/3 + f₂/2 (acoustic contrast factor)
//! - f₁ = 1 − κ_p/κ_m (compressibility contrast)
//! - f₂ = 2(ρ_p − ρ_m)/(2ρ_p + ρ_m) (density contrast)
//!
//! For blood cells in plasma (ρ_m = 1025 kg/m³, κ_m = 4.5e-10 Pa⁻¹):
//! - CTCs (ρ ≈ 1050 kg/m³, κ ≈ 4.2e-10 Pa⁻¹): Φ ≈ 0.030 (positive → moves to node)
//! - RBCs (ρ ≈ 1090 kg/m³, κ ≈ 3.4e-10 Pa⁻¹): Φ ≈ 0.102 (positive → moves to node)
//! - WBCs (ρ ≈ 1060 kg/m³, κ ≈ 3.9e-10 Pa⁻¹): Φ ≈ 0.056
//!
//! **Reference**: Gor'kov, L.P. (1962). "On the forces acting on a small
//! particle in an acoustical field in an ideal fluid",
//! *Sov. Phys. Dokl.* 6:773-775.
//! Also: Bruus, H. (2012). "Acoustofluidics 7: The acoustic radiation force
//! on small particles", *Lab Chip* 12:1014-1021.

use std::f64::consts::PI;

// ── Typical blood-cell properties in plasma ────────────────────────────────

/// Typical CTC density \[kg/m³\].
pub const RHO_CTC: f64 = 1050.0;
/// Typical CTC compressibility \[Pa⁻¹\].
pub const KAPPA_CTC: f64 = 4.2e-10;

/// Typical RBC density \[kg/m³\].
pub const RHO_RBC: f64 = 1090.0;
/// Typical RBC compressibility \[Pa⁻¹\].
pub const KAPPA_RBC: f64 = 3.4e-10;

/// Typical WBC density \[kg/m³\].
pub const RHO_WBC: f64 = 1060.0;
/// Typical WBC compressibility \[Pa⁻¹\].
pub const KAPPA_WBC: f64 = 3.9e-10;

/// Blood plasma density \[kg/m³\].
pub const RHO_PLASMA: f64 = 1025.0;
/// Blood plasma compressibility \[Pa⁻¹\].
pub const KAPPA_PLASMA: f64 = 4.5e-10;
/// Speed of sound in blood plasma \[m/s\].
pub const SPEED_OF_SOUND_PLASMA: f64 = 1540.0;

// ── Public API ─────────────────────────────────────────────────────────────

/// Acoustic contrast factor Φ for a sphere in a fluid.
///
/// The contrast factor determines whether a particle migrates toward pressure
/// nodes (Φ > 0) or antinodes (Φ < 0) in a standing acoustic wave.
///
/// ```text
/// Φ = f₁/3 + f₂/2
/// f₁ = 1 − κ_p/κ_m   (compressibility contrast)
/// f₂ = 2(ρ_p − ρ_m)/(2ρ_p + ρ_m)   (density contrast)
/// ```
///
/// # Arguments
///
/// * `rho_particle` — particle density \[kg/m³\]
/// * `rho_medium` — medium density \[kg/m³\]
/// * `kappa_particle` — particle compressibility \[Pa⁻¹\]
/// * `kappa_medium` — medium compressibility \[Pa⁻¹\]
///
/// # Example
///
/// ```
/// use cfd_1d::physics::hemolysis::acoustic_radiation::*;
/// let phi = acoustic_contrast_factor(RHO_CTC, RHO_PLASMA, KAPPA_CTC, KAPPA_PLASMA);
/// assert!(phi > 0.0, "CTCs should migrate toward pressure nodes");
/// ```
#[inline]
#[must_use]
pub fn acoustic_contrast_factor(
    rho_particle: f64,
    rho_medium: f64,
    kappa_particle: f64,
    kappa_medium: f64,
) -> f64 {
    let f1 = 1.0 - kappa_particle / kappa_medium;
    let f2 = 2.0 * (rho_particle - rho_medium) / (2.0 * rho_particle + rho_medium);
    f1 / 3.0 + f2 / 2.0
}

/// Acoustic radiation force magnitude on a sphere in a standing wave \[N\].
///
/// ```text
/// F_rad = (4π/3) a³ · Φ · k · E_ac · sin(2kx)
/// ```
///
/// The force is position-dependent: zero at pressure nodes (x = 0) and maximal
/// midway between nodes (x = λ/8).
///
/// # Arguments
///
/// * `radius_m` — cell radius \[m\]
/// * `contrast_factor` — Φ (from [`acoustic_contrast_factor`])
/// * `frequency_hz` — ultrasound frequency \[Hz\]
/// * `pressure_amplitude_pa` — acoustic pressure amplitude p₀ \[Pa\]
/// * `rho_medium` — medium density \[kg/m³\]
/// * `speed_of_sound` — speed of sound in medium \[m/s\]
/// * `position_x` — position relative to nearest pressure node \[m\]
///
/// # Example
///
/// ```
/// use cfd_1d::physics::hemolysis::acoustic_radiation::*;
/// let phi = acoustic_contrast_factor(RHO_RBC, RHO_PLASMA, KAPPA_RBC, KAPPA_PLASMA);
/// let f = acoustic_radiation_force(3.5e-6, phi, 2.0e6, 500_000.0, RHO_PLASMA, SPEED_OF_SOUND_PLASMA, 0.0);
/// assert!(f.abs() < 1e-30, "Force should be zero at a pressure node");
/// ```
#[inline]
#[must_use]
pub fn acoustic_radiation_force(
    radius_m: f64,
    contrast_factor: f64,
    frequency_hz: f64,
    pressure_amplitude_pa: f64,
    rho_medium: f64,
    speed_of_sound: f64,
    position_x: f64,
) -> f64 {
    let k = 2.0 * PI * frequency_hz / speed_of_sound;
    let e_ac = pressure_amplitude_pa.powi(2) / (4.0 * rho_medium * speed_of_sound.powi(2));
    let volume = (4.0 / 3.0) * PI * radius_m.powi(3);
    volume * contrast_factor * k * e_ac * (2.0 * k * position_x).sin()
}

/// Acoustic energy density \[J/m³\].
///
/// ```text
/// E_ac = p₀² / (4ρc²)
/// ```
///
/// # Arguments
///
/// * `pressure_amplitude_pa` — acoustic pressure amplitude p₀ \[Pa\]
/// * `rho_medium` — medium density \[kg/m³\]
/// * `speed_of_sound` — speed of sound \[m/s\]
#[inline]
#[must_use]
pub fn acoustic_energy_density(
    pressure_amplitude_pa: f64,
    rho_medium: f64,
    speed_of_sound: f64,
) -> f64 {
    pressure_amplitude_pa.powi(2) / (4.0 * rho_medium * speed_of_sound.powi(2))
}

// ── Tests ──────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    // ── acoustic_contrast_factor ───────────────────────────────────────────

    #[test]
    fn test_contrast_factor_ctc_positive() {
        let phi = acoustic_contrast_factor(RHO_CTC, RHO_PLASMA, KAPPA_CTC, KAPPA_PLASMA);
        // CTC in plasma: Φ ≈ 0.030 (positive → migrates to pressure node)
        assert!(phi > 0.0, "CTC contrast factor must be positive, got {phi}");
        assert!(
            (phi - 0.030).abs() < 0.01,
            "CTC contrast factor should be approximately 0.030, got {phi}"
        );
    }

    #[test]
    fn test_contrast_factor_rbc() {
        let phi = acoustic_contrast_factor(RHO_RBC, RHO_PLASMA, KAPPA_RBC, KAPPA_PLASMA);
        // RBC in plasma: Φ ≈ 0.102 (positive → migrates to pressure node)
        assert!(phi > 0.0, "RBC contrast factor must be positive, got {phi}");
        assert!(
            (phi - 0.102).abs() < 0.01,
            "RBC contrast factor should be approximately 0.102, got {phi}"
        );
    }

    #[test]
    fn test_contrast_factor_wbc() {
        let phi = acoustic_contrast_factor(RHO_WBC, RHO_PLASMA, KAPPA_WBC, KAPPA_PLASMA);
        // WBC in plasma: Φ ≈ 0.056
        assert!(phi > 0.0, "WBC contrast factor must be positive, got {phi}");
        assert!(
            (phi - 0.056).abs() < 0.01,
            "WBC contrast factor should be approximately 0.056, got {phi}"
        );
    }

    // ── acoustic_radiation_force ───────────────────────────────────────────

    #[test]
    fn test_radiation_force_zero_at_node() {
        let phi = acoustic_contrast_factor(RHO_RBC, RHO_PLASMA, KAPPA_RBC, KAPPA_PLASMA);
        let f = acoustic_radiation_force(
            3.5e-6, // RBC radius ~3.5 µm
            phi,
            2.0e6,     // 2 MHz
            500_000.0, // 500 kPa
            RHO_PLASMA,
            SPEED_OF_SOUND_PLASMA,
            0.0, // at pressure node
        );
        assert!(
            f.abs() < 1e-30,
            "Force at pressure node should be zero, got {f}"
        );
    }

    #[test]
    fn test_radiation_force_max_between_nodes() {
        let phi = acoustic_contrast_factor(RHO_RBC, RHO_PLASMA, KAPPA_RBC, KAPPA_PLASMA);
        let freq = 2.0e6_f64;
        let lambda = SPEED_OF_SOUND_PLASMA / freq;
        // At x = λ/8, sin(2kx) = sin(2 · 2π/λ · λ/8) = sin(π/2) = 1 → maximum force
        let x_max = lambda / 8.0;
        let f_max = acoustic_radiation_force(
            3.5e-6,
            phi,
            freq,
            500_000.0,
            RHO_PLASMA,
            SPEED_OF_SOUND_PLASMA,
            x_max,
        );
        assert!(
            f_max > 0.0,
            "Force at λ/8 should be positive (max), got {f_max}"
        );

        // Verify it is indeed the maximum by checking nearby positions
        let f_nearby = acoustic_radiation_force(
            3.5e-6,
            phi,
            freq,
            500_000.0,
            RHO_PLASMA,
            SPEED_OF_SOUND_PLASMA,
            x_max * 0.5,
        );
        assert!(
            f_max > f_nearby,
            "Force at λ/8 ({f_max}) should exceed force at λ/16 ({f_nearby})"
        );
    }

    #[test]
    fn test_radiation_force_scales_with_radius_cubed() {
        let phi = acoustic_contrast_factor(RHO_RBC, RHO_PLASMA, KAPPA_RBC, KAPPA_PLASMA);
        let freq = 2.0e6;
        let lambda = SPEED_OF_SOUND_PLASMA / freq;
        let x = lambda / 8.0;
        let r1 = 3.5e-6;
        let r2 = 2.0 * r1; // double the radius

        let f1 = acoustic_radiation_force(
            r1,
            phi,
            freq,
            500_000.0,
            RHO_PLASMA,
            SPEED_OF_SOUND_PLASMA,
            x,
        );
        let f2 = acoustic_radiation_force(
            r2,
            phi,
            freq,
            500_000.0,
            RHO_PLASMA,
            SPEED_OF_SOUND_PLASMA,
            x,
        );

        let ratio = f2 / f1;
        assert!(
            (ratio - 8.0).abs() < 1e-10,
            "Doubling radius should give 8x force (r³ scaling), got ratio {ratio}"
        );
    }

    // ── acoustic_energy_density ────────────────────────────────────────────

    #[test]
    fn test_energy_density_positive() {
        let e = acoustic_energy_density(500_000.0, RHO_PLASMA, SPEED_OF_SOUND_PLASMA);
        assert!(
            e > 0.0,
            "Energy density must be positive for p₀ > 0, got {e}"
        );
    }

    #[test]
    fn test_energy_density_scales_with_pressure_squared() {
        let e1 = acoustic_energy_density(100_000.0, RHO_PLASMA, SPEED_OF_SOUND_PLASMA);
        let e2 = acoustic_energy_density(200_000.0, RHO_PLASMA, SPEED_OF_SOUND_PLASMA);
        let ratio = e2 / e1;
        assert!(
            (ratio - 4.0).abs() < 1e-10,
            "Doubling pressure should give 4x energy density, got ratio {ratio}"
        );
    }

    #[test]
    fn test_energy_density_zero_pressure() {
        let e = acoustic_energy_density(0.0, RHO_PLASMA, SPEED_OF_SOUND_PLASMA);
        assert!(
            e.abs() < 1e-30,
            "Energy density should be zero for zero pressure, got {e}"
        );
    }
}
